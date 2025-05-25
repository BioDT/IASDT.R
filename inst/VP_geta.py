# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# - Authors: Ahmed El-Gabbas [elgabbas[at]outlook[dot]com]
# - Last updated: 24.05.2025
# - Affiliation: Helmholtz Centre for Environmental Research - UFZ, Germany
# - License: MIT License
#
# - Description:
#       
#   This Python script helps in computing variance partitioning for
#   Hierarchical Modelling of Species Communities (Hmsc) models. It serves as a 
#   backend utility called within the R function 
#   `IASDT.R::variance_partitioning_compute` from the `IASDT.R` package, which 
#   models the distribution of invasive alien plant species in Europe. The 
#   script adapts the functionality of the R function 
#   Hmsc::computeVariancePartitioning`, leveraging TensorFlow for efficient 
#   matrix multiplications. The script includes robust error handling, memory 
#   management, and file validation to ensure reliability and efficiency.
#   
#   This script replaces the `geta` function defined inside
#   the `computeVariancePartitioning` function in the R package `Hmsc`:
#
#   geta = function(a){
#     switch(
#       class(hM$X)[1L],
#       matrix = {res = hM$X%*%(t(hM$Tr%*%t(a$Gamma)))},
#       list = {
#         res = matrix(NA,hM$ny,hM$ns)
#         for(j in 1:hM$ns) res[,j] =  hM$X[[j]]%*%(t(hM$Tr[j,]%*%t(a$Gamma)))
#       }
#     )
#     return(res)
#   }
#
# - CITATION: El-Gabbas, A. (2025). IASDT.R: Modelling the distribution of invasive alien plant species in Europe. https://doi.org/10.5281/zenodo.1483438, R package version 0.1.X; https://biodt.github.io/IASDT.R/.
# 
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


import multiprocessing
import argparse
import logging
import warnings
import traceback
import os
import numpy as np
import pandas as pd
import pyarrow.feather as feather
import tensorflow as tf
import gc

# ======================================================================
# Disable Additional Logs and Configure TensorFlow
# ======================================================================

# Suppress TensorFlow informational and warning messages, showing only errors
# https://stackoverflow.com/questions/55081911/
logging.getLogger('tensorflow').setLevel(logging.ERROR)

# Set TensorFlow logging level to show only errors
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Disable oneDNN optimizations to reduce unnecessary logs
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

# Suppress all TensorFlow warnings
warnings.filterwarnings("ignore", category=UserWarning, module="tensorflow")

# Clear any existing TensorFlow session to start fresh
tf.keras.backend.clear_session()
# Set TensorFlow logger to display only error-level messages
tf.get_logger().setLevel('ERROR')

# https://www.tensorflow.org/api_docs/python/tf/compat/v1/reset_default_graph
# tf.compat.v1.reset_default_graph()
# tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

# ======================================================================
# TensorFlow and Environment Configuration
# ======================================================================

# Detect available GPUs for TensorFlow computations
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        # Configure each GPU to allocate memory dynamically, preventing full 
        # reservation
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        # Set a fixed memory limit (e.g., 64 GB) to avoid exceeding the GPU memory
        #tf.config.experimental.set_virtual_device_configuration(
            #gpus[0],
            #[tf.config.experimental.VirtualDeviceConfiguration(memory_limit=64000)]
        #)
    except RuntimeError as e:
        print(f"Error setting memory configurations: {e}")

# ======================================================================
# Function Definitions
# ======================================================================

def load_feather_tensor(file_path, dtype):
    """
    Load a Feather file and return its data as a TensorFlow tensor.

    This function reads a Feather file into a NumPy array, ensures memory 
    contiguity, converts it to a TensorFlow tensor, and manages memory by 
    cleaning up intermediate objects.

    Args:
        file_path (str): Path to the Feather file.
        dtype (tf.DType): TensorFlow data type (e.g., tf.float32 or tf.float64).

    Returns:
        tf.Tensor: Data loaded as a TensorFlow tensor.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If Feather data reading fails.
    """
    
    # Verify that the input file exists before attempting to load
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: The specified file '{file_path}' does not exist.")
    try:
        # Load Feather file data into a NumPy array via pandas
        feather_data = feather.read_feather(file_path).to_numpy()
        # Ensure the array is contiguous in memory for optimal performance
        feather_data = np.ascontiguousarray(feather_data)
        # Convert the NumPy array to a TensorFlow tensor with specified data type
        tensor = tf.convert_to_tensor(feather_data, dtype=dtype)
        # Delete the intermediate NumPy array to free memory
        del feather_data
        # Trigger garbage collection to reclaim memory immediately
        gc.collect()
        return tensor
    except ValueError as e:
        raise ValueError(f"Failed to read Feather file at {file_path}: {e}")

# ======================================================================
# ======================================================================

def process_column(tr, x, gamma_chunk, col_offset, file_output):
    """
    Process a column chunk of the gamma matrix and compute x %*% (tr %*% gamma).

    This function iterates over columns in a gamma chunk, performs matrix 
    multiplications, saves results to Feather files, and validates them. It 
    skips existing valid files and retries failed operations up to 5 times.

    Args:
        tr (tf.Tensor): Transformation matrix.
        x (tf.Tensor): Data matrix.
        gamma_chunk (tf.Tensor): Chunk of the gamma matrix to process.
        col_offset (int): Column offset for global indexing.
        file_output (str): Base path for output Feather files.
    """
    try:
        # Iterate over each column in the gamma chunk
        for col in range(gamma_chunk.shape[1]):
            # Compute the global column index based on the offset
            col_index = col + col_offset
            # Format the file index as a 4-digit padded string
            file_index = f"{col_index + 1:04d}"
            # Construct the output Feather file path with the index
            feather_file = f"{file_output.replace('.feather', f'_{file_index}.feather')}"
            # Check if the output file already exists
            if os.path.exists(feather_file):
                try:
                    # Attempt to read the file to verify its validity
                    feather.read_feather(feather_file)
                    # Skip processing if the file is valid
                    continue
                except Exception:
                    # Remove the file if it exists but is corrupted
                    os.remove(feather_file)
            # Define maximum retries for handling transient errors
            max_retries = 5
            # Attempt processing with retries
            for attempt in range(max_retries):
                try:
                    # Extract and reshape the current gamma column into a 
                    # column vector
                    tr_gamma = tf.reshape(gamma_chunk[:, col], [-1, 1])
                    # Compute the matrix product tr %*% gamma (transposed for 
                    # correct orientation)
                    tr_gamma = tf.matmul(tr, tr_gamma, transpose_b=True)
                    # Compute the matrix product x %*% (tr %*% gamma)
                    res = tf.matmul(x, tr_gamma, transpose_b=True)
                    # Convert the TensorFlow result to a pandas DataFrame
                    res_df = pd.DataFrame(res.numpy())
                    # Save the result to a Feather file
                    res_df.to_feather(feather_file)
                    # Validate the saved file by reading it back
                    feather.read_feather(feather_file)
                    # Clean up temporary variables to free memory
                    del tr_gamma, res, res_df
                    # Trigger garbage collection after successful save
                    gc.collect()
                    break  # Exit retry loop on success
                except Exception as e:
                    # Remove any corrupted file from this attempt
                    if os.path.exists(feather_file):
                        os.remove(feather_file)
                    # Check if this was the final attempt
                    if attempt == max_retries - 1:
                        print(f"Error: Failed to save valid Feather file {feather_file} after {max_retries} attempts: {e}")
                        raise
                    # Clean up temporary variables even on failure
                    del tr_gamma, res, res_df
                    # Trigger garbage collection before retrying
                    gc.collect()
    except Exception as e:
        print(f"Job failed with exception: {e}")
        raise

# ======================================================================
# ======================================================================

def chunkify(data, chunk_size):
    """
    Split data into chunks of specified size.

    This generator function yields chunks of the input tensor along with their 
    starting indices, facilitating sequential processing.

    Args:
        data (tf.Tensor): Input data tensor to be chunked.
        chunk_size (int): Number of columns per chunk.

    Yields:
        tuple: A chunk of the data (tf.Tensor) and its starting index (int).
    """
    # Iterate over the column dimension in steps of chunk_size
    for start in range(0, data.shape[1], chunk_size):
        # Yield a slice of the data and its starting index
        yield data[:, start:start + chunk_size], start

# ======================================================================
# ======================================================================

def geta(file_tr, file_x, file_gamma, use_single, file_output, chunk_size):
    """
    Perform sequential matrix multiplication and save results.

    This function loads input matrices, splits the gamma matrix into chunks, 
    and processes each chunk sequentially, saving results to Feather files.

    Args:
        file_tr (str): Path to the tr matrix Feather file.
        file_x (str): Path to the x matrix Feather file.
        file_gamma (str): Path to the gamma matrix Feather file.
        use_single (bool): If True, use single precision (float32); otherwise, double (float64).
        file_output (str): Base path for output Feather files.
        chunk_size (int): Number of columns to process per chunk.
    """
    # Validate that chunk_size is positive
    if chunk_size <= 0:
        raise ValueError("Chunk size must be a positive integer.")
    
    # Select data type based on precision preference
    dtype = tf.float32 if use_single else tf.float64
    # Load the transformation matrix from its Feather file
    tr = load_feather_tensor(file_tr, dtype)
    # Load the data matrix from its Feather file
    x = load_feather_tensor(file_x, dtype)
    # Load the gamma matrix from its Feather file
    gamma = load_feather_tensor(file_gamma, dtype)
    # Free memory after loading large tensors
    gc.collect()
    
    # Create a list of tasks by splitting gamma into chunks
    tasks = [(tr, x, chunk, offset, file_output) for chunk, offset in chunkify(gamma, chunk_size)]
    # Delete the gamma tensor to free memory since chunks are now in tasks
    del gamma
    # Trigger garbage collection after deleting gamma
    gc.collect()
    
    # Process each chunk sequentially
    for args in tasks:
        try:
            # Process the current chunk and save results
            process_column(*args)
            # Delete task arguments to free memory
            del args
            # Trigger garbage collection after processing
            gc.collect()
        except Exception as e:
            print(f"Error processing task: {e}")
            traceback.print_exc()

# ======================================================================
# Main Function
# ======================================================================

def main():
    """
    Parse command-line arguments and initiate the computation.

    This function sets up an argument parser, collects user inputs, and calls 
    the `geta` function to perform the matrix operations.
    """
    # Initialize the argument parser with a description
    parser = argparse.ArgumentParser(description="Sequential TensorFlow-based matrix operations.")
    # Add required argument for the tr matrix file path
    parser.add_argument("--tr", required=True, help="Path to the tr Feather file")
    # Add required argument for the x matrix file path
    parser.add_argument("--x", required=True, help="Path to the x Feather file")
    # Add required argument for the gamma matrix file path
    parser.add_argument("--gamma", required=True, help="Path to the gamma Feather file")
    # Add required argument for the output file path
    parser.add_argument("--output", required=True, help="Output Feather file path")
    # Add optional flag for single precision computation
    parser.add_argument("--use_single", action="store_true", help="Use single precision for computations.")
    # Add optional argument for chunk size with a default value
    parser.add_argument("--chunk_size", type=int, default=20, help="Size of chunks for processing")
    # Parse the command-line arguments
    args = parser.parse_args()
    # Execute the matrix computation with provided arguments
    geta(args.tr, args.x, args.gamma, args.use_single, args.output, args.chunk_size)

if __name__ == "__main__":
    # Set multiprocessing start method to 'spawn' for compatibility across platforms
    multiprocessing.set_start_method("spawn", force=True)
    print("Done")
    # Run the main function to start the script
    main()
