# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# - Authors: Ahmed El-Gabbas [elgabbas[at]outlook[dot]com]
# - Last updated: 24.05.2025
# - Affiliation: Helmholtz Centre for Environmental Research - UFZ, Germany
# - License: MIT License
#
# - Description:
#   
#   This Python script computes variance partitioning for Hierarchical 
#   Modelling of Species Communities (Hmsc) models. It replaces the R function 
#   `gemu` (see below) from the `Hmsc::computeVariancePartitioning` R function 
#   and is designed to be called within 
#   `IASDT.R::variance_partitioning_compute` from the `IASDT.R` package, which 
#   models the distribution of invasive alien plant species in Europe. The 
#   script uses TensorFlow for efficient matrix multiplications and processes 
#   data in parallel with the `loky` library, handling large datasets by 
#   splitting the Gamma matrix into chunks. Results are saved as individual 
#   Feather files. It includes error handling and memory management for 
#   stability.
#
#
#   gemu = function(a){
#     res = t(hM$Tr%*%t(a$Gamma))
#     return(res)
#   }
# 
# - CITATION: El-Gabbas, A. (2025). IASDT.R: Modelling the distribution of 
#   invasive alien plant species in Europe. https://doi.org/10.5281/zenodo.1483438, 
#   R package version 0.1.X; https://biodt.github.io/IASDT.R/.
# 
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import multiprocessing
import argparse
import logging
import warnings
import os
import traceback
import numpy as np
import pandas as pd
import pyarrow.feather as feather
import tensorflow as tf
import gc
from loky import get_reusable_executor

# ======================================================================
# Disable Additional Logs and Configure TensorFlow
# ======================================================================

# Suppress TensorFlow logs, showing only errors
# https://stackoverflow.com/questions/55081911/
logging.getLogger('tensorflow').setLevel(logging.ERROR)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# Disable oneDNN optimizations to reduce log noise
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
# Ignore TensorFlow user warnings
warnings.filterwarnings("ignore", category=UserWarning, module="tensorflow")
# Clear TensorFlow session for a fresh start
tf.keras.backend.clear_session()
# Set TensorFlow logger to error level only
# https://www.tensorflow.org/api_docs/python/tf/get_logger
tf.get_logger().setLevel('ERROR')

# https://www.tensorflow.org/api_docs/python/tf/compat/v1/reset_default_graph
# tf.compat.v1.reset_default_graph()
# tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

# ======================================================================
# TensorFlow and Environment Configuration
# ======================================================================

# Detect available GPUs
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        # Enable dynamic memory growth for GPUs
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
    Load a Feather file into a TensorFlow tensor.
    
    Reads data into a NumPy array, ensures contiguous memory, and converts it 
    to a TensorFlow tensor with the specified data type.
    
    Args:
        file_path (str): Path to the Feather file.
        dtype (tf.DType): TensorFlow data type (e.g., tf.float32 or tf.float64).
    
    Returns:
        tf.Tensor: Loaded data as a TensorFlow tensor.
    
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If reading the Feather file fails.
    """

    # Check if the file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    try:
        # Load Feather data into a NumPy array
        data = feather.read_feather(file_path).to_numpy()
        # Ensure contiguous memory layout
        data = np.ascontiguousarray(data)
        # Convert to TensorFlow tensor
        tensor = tf.convert_to_tensor(data, dtype=dtype)
        # Free memory by deleting NumPy array
        del data
        gc.collect()
        return tensor
    except ValueError as e:
        raise ValueError(f"Error reading {file_path}: {e}")

# ======================================================================
# ======================================================================

def process_column(tr, gamma_chunk, col_offset, file_output):
    """
    Process a chunk of the Gamma matrix and save results.
    
    Computes tr %*% t(gamma) for each column in the chunk and saves the result 
    as a Feather file. Skips processing if the output file already exists.
    
    Args:
        tr (tf.Tensor): Transformation matrix.
        gamma_chunk (tf.Tensor): Chunk of the Gamma matrix.
        col_offset (int): Offset for global column indexing.
        file_output (str): Base path for output Feather files.
    
    Raises:
        Exception: If processing fails.
    """

    try:
        # Iterate over each column in the chunk
        for col in range(gamma_chunk.shape[1]):
            # Calculate global column index
            col_index = col + col_offset
            # Generate padded file index (e.g., 0001)
            file_index = f"{col_index + 1:04d}"
            # Construct output file path
            feather_file = file_output.replace('.feather', f'_{file_index}.feather')
            # Skip if output file already exists
            if os.path.exists(feather_file):
                continue
            # Extract and reshape the column to a column vector
            gamma_col = tf.reshape(gamma_chunk[:, col], [-1, 1])
            # Compute tr %*% t(gamma_col)
            tr_gamma = tf.matmul(tr, gamma_col, transpose_b=True)
            # Transpose and convert to DataFrame
            res = pd.DataFrame(tf.transpose(tr_gamma).numpy())
            # Save to Feather file
            res.to_feather(feather_file)
            # Clean up memory
            del gamma_col, tr_gamma, res
            gc.collect()
    except Exception as e:
        print(f"Processing failed: {e}")
        raise

# ======================================================================
# ======================================================================

def chunkify(data, chunk_size):
    """
    Split a tensor into chunks along the column dimension.
    
    Yields chunks of the specified size and their starting indices.
    
    Args:
        data (tf.Tensor): Input tensor to split.
        chunk_size (int): Number of columns per chunk.
    
    Yields:
        tuple: (chunk, start_index)
    """

    # Iterate over columns in steps of chunk_size
    for start in range(0, data.shape[1], chunk_size):
        # Extract chunk and yield with start index
        yield data[:, start:start + chunk_size], start

# ======================================================================
# ======================================================================

def gemu(file_tr, file_gamma, use_single, file_output, ncores, chunk_size):
    """
    Perform parallel matrix operations and save results.
    
    Loads tr and Gamma matrices, splits Gamma into chunks, and processes each 
    chunk in parallel using loky. Results are saved as Feather files.
    
    Args:
        file_tr (str): Path to tr matrix Feather file.
        file_gamma (str): Path to Gamma matrix Feather file.
        use_single (bool): Use float32 if True, else float64.
        file_output (str): Base path for output Feather files.
        ncores (int): Number of cores for parallel processing.
        chunk_size (int): Columns per chunk.
    """
    # Validate chunk_size and ncores
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if ncores <= 0:
        raise ValueError("ncores must be positive")
    
    # Set data type based on precision flag
    dtype = tf.float32 if use_single else tf.float64
    # Load tr and Gamma matrices
    tr = load_feather_tensor(file_tr, dtype)
    gamma = load_feather_tensor(file_gamma, dtype)
    
    # Create tasks by splitting Gamma into chunks
    tasks = [(tr, chunk, offset, file_output) for chunk, offset in chunkify(gamma, chunk_size)]
    # Free memory by deleting Gamma
    del gamma
    gc.collect()
    
    # Process tasks in parallel with loky
    with get_reusable_executor(max_workers=ncores, timeout=600) as executor:
        # Submit tasks to executor
        futures = [executor.submit(process_column, *args) for args in tasks]
        # Check results and handle errors
        for future in futures:
            try:
                future.result()
            except Exception as e:
                print(f"Worker error: {e}")
                traceback.print_exc()

# ======================================================================
# Main Function
# ======================================================================

def main():
    """
    Parse arguments and start the computation.
    
    Sets up command-line arguments and initiates the matrix operations.
    """

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Parallel TensorFlow matrix operations")
    parser.add_argument("--tr", required=True, help="Path to tr Feather file")
    parser.add_argument("--gamma", required=True, help="Path to Gamma Feather file")
    parser.add_argument("--output", required=True, help="Output file base path")
    parser.add_argument("--use_single", action="store_true", help="Use float32")
    parser.add_argument("--ncores", type=int, default=4, help="Number of cores")
    parser.add_argument("--chunk_size", type=int, default=20, help="Chunk size")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Start computation
    gemu(args.tr, args.gamma, args.use_single, args.output, args.ncores, args.chunk_size)

# ======================================================================
# ======================================================================

if __name__ == "__main__":
    # Set spawn method for cross-platform compatibility
    multiprocessing.set_start_method("spawn", force=True)
    print("Done")
    # Execute main function
    main()
