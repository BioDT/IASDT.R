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
#   Hierarchical Modelling of Species Communities (Hmsc) models. It replaces 
#   the R function `getf` (see below) from the 
#   `Hmsc::computeVariancePartitioning` R function and is designed to be called 
#   within `IASDT.R::variance_partitioning_compute` from the `IASDT.R` package, 
#   which models the distribution of invasive alien plant species in Europe.
#   The script uses TensorFlow for efficient matrix multiplications and 
#   processes data in parallel with the `loky` library, handling large datasets 
#   by splitting the Beta matrix into chunks stored as Feather files. Results 
#   are saved as individual Feather files. It includes robust error handling 
#   and memory management for reliability.
#
#   getf = function(a){
#     switch(
#       class(hM$X)[1L],
#       matrix = {res = hM$X%*%(a$Beta)},
#       list = {
#         res = matrix(NA,hM$ny,hM$ns)
#         for(j in 1:hM$ns) res[,j] = hM$X[[j]]%*%a$Beta[,j]
#       }
#     )
#     return(res)
#   }

# - CITATION: El-Gabbas, A. (2025). IASDT.R: Modelling the distribution of 
#   invasive alien plant species in Europe. https://doi.org/10.5281/zenodo.14834384, 
#   R package version 0.1.X; https://biodt.github.io/IASDT.R/.
# 
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import multiprocessing
import argparse
import logging
import warnings
import os
import re
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
# Set TensorFlow logging level to show only errors
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

def configure_gpu_memory(ncores, total_memory_mb=64 * 1024):
    """
    Configure GPU memory limits based on the number of cores.
    
    Allocates memory per worker to prevent exceeding total GPU memory. The 
    total memory (default 64GB) is divided equally among the specified cores.
    
    Args:
        ncores (int): Number of cores for parallel processing.
        total_memory_mb (int): Total GPU memory in MB (default: 64GB).
    
    Raises:
        ValueError: If ncores is not positive.
        RuntimeError: If memory configuration fails.
    """

    # Get list of available GPUs
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                # Ensure ncores is valid
                if ncores <= 0:
                    raise ValueError("ncores must be greater than zero.")
                # Calculate memory per worker
                memory_limit = total_memory_mb // ncores
                
                # Set memory limit for the GPU
                tf.config.set_logical_device_configuration(
                    gpu,
                    [tf.config.LogicalDeviceConfiguration(memory_limit=memory_limit)]
                )
                # Enable memory growth
                # tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(f"Error setting memory limits: {e}")

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

def process_beta_file(x, beta_file, task_index, dtype, file_output):
    """
    Perform matrix multiplication for a Beta chunk and save the result.
    
    Loads a Beta chunk, computes x %*% Beta, and saves the result as a Feather 
    file. Includes retry logic (up to 5 attempts) for robustness.
    
    Args:
        x (tf.Tensor): Input matrix x.
        beta_file (str): Path to the Beta chunk Feather file.
        task_index (int): Index for naming the output file.
        dtype (tf.DType): Data type for computations.
        file_output (str): Base path for the output Feather file.
    
    Raises:
        ValueError: If matrix dimensions do not match.
        Exception: If processing fails after retries.
    """

    try:
        # Generate output file name with padded index
        file_index = f"{task_index + 1:04d}"
        feather_file = file_output.replace(".feather", f"_{file_index}.feather")
        
        # Skip if output file already exists
        if os.path.exists(feather_file):
            return
        
        # Load Beta chunk from Feather file
        beta_chunk = load_feather_tensor(beta_file, dtype)
        
        # Validate matrix dimensions
        if beta_chunk.shape[0] != x.shape[1]:
            raise ValueError(f"x cols ({x.shape[1]}) != Beta rows ({beta_chunk.shape[0]})")
        
        # Retry up to 5 times on failure
        max_retries = 5
        for attempt in range(1, max_retries + 1):
            try:
                # Compute matrix multiplication
                res = tf.matmul(x, beta_chunk)
                # Convert result to NumPy array
                result_array = res.numpy()
                # Create DataFrame with string column names
                num_cols = result_array.shape[1]
                df = pd.DataFrame(result_array, columns=[str(i) for i in range(num_cols)])
                # Save to Feather file
                df.to_feather(feather_file)
                # Clean up memory
                del result_array, df
                gc.collect()
                return  # Success, exit loop
            except Exception as e:
                if attempt == max_retries:
                    print(f"Max retries reached for {feather_file}")
                    raise
                # Clean up before retry
                if 'result_array' in locals():
                    del result_array
                if 'df' in locals():
                    del df
                gc.collect()
    except Exception as e:
        print(f"Processing failed: {e}")
        raise

# ======================================================================
# ======================================================================

def validate_and_get_beta_files(beta_dir, prefix="vp_beta_", extension=".feather"):
    """
    Validate and list Beta files in the directory.
    
    Ensures all files matching '{prefix}{index}{extension}' exist with 
    continuous indices. Returns a sorted list of file paths.
    
    Args:
        beta_dir (str): Directory with Beta files.
        prefix (str): File name prefix (default: "vp_beta_").
        extension (str): File extension (default: ".feather").
    
    Returns:
        list: Sorted list of Beta file paths.
    
    Raises:
        FileNotFoundError: If files are missing or none are found.
    """

    # Define pattern for Beta files (e.g., vp_beta_0001.feather)
    pattern = re.compile(rf"^{re.escape(prefix)}(\d+){re.escape(extension)}$")
    beta_files = []
    file_indices = []
    
    # Scan directory for matching files
    for f in os.listdir(beta_dir):
        match = pattern.match(f)
        if match:
            beta_files.append(os.path.join(beta_dir, f))
            file_indices.append(int(match.group(1)))
    
    # Check if any files were found
    if not beta_files:
        raise FileNotFoundError(f"No {prefix}*{extension} files in {beta_dir}")
    
    # Sort indices and check for gaps
    file_indices = sorted(file_indices)
    expected = list(range(file_indices[0], file_indices[-1] + 1))
    missing = set(expected) - set(file_indices)
    
    # Report missing files if any
    if missing:
        missing_str = ", ".join([f"{prefix}{i:04d}{extension}" for i in sorted(missing)])
        raise FileNotFoundError(f"Missing files: {missing_str}")
    
    # Return sorted file list
    return sorted(beta_files)

# ======================================================================
# ======================================================================

def getf(file_x, beta_dir, use_single, file_output, ncores):
    """
    Perform matrix multiplication on Beta files and save results.
    
    Loads the x matrix, processes Beta files either sequentially (ncores=1) or 
    in parallel (ncores>1), and saves results as Feather files.
    
    Args:
        file_x (str): Path to the x matrix Feather file.
        beta_dir (str): Directory with Beta Feather files.
        use_single (bool): Use float32 if True, else float64.
        file_output (str): Base path for output Feather files.
        ncores (int): Number of cores for parallel processing.
    """

    # Set data type based on precision flag
    dtype = tf.float32 if use_single else tf.float64
    # Load x matrix from Feather file
    x = load_feather_tensor(file_x, dtype)
    
    # Get and validate Beta files
    beta_files = validate_and_get_beta_files(beta_dir)
    
    try:
        if ncores == 1:
            # Process files sequentially
            for i, beta_file in enumerate(beta_files):
                process_beta_file(x, beta_file, i, dtype, file_output)
        else:
            # Process files in parallel with loky
            with get_reusable_executor(max_workers=ncores, timeout=600) as executor:
                # Submit tasks for parallel execution
                futures = [
                    executor.submit(process_beta_file, x, beta_file, i, dtype, file_output)
                    for i, beta_file in enumerate(beta_files)
                ]
                # Check results and handle errors
                for future in futures:
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Worker error: {e}")
                        traceback.print_exc()
    except Exception as e:
        print(f"Parallel processing error: {e}")
        raise

# ======================================================================
# Main Function
# ======================================================================

def main():
    """
    Parse arguments and start the computation.
    
    Sets up command-line arguments, configures GPU memory, and initiates the 
    matrix multiplication process.
    """

    # Set up argument parser
    parser = argparse.ArgumentParser(description="TensorFlow matrix operations")
    parser.add_argument("--x", required=True, help="Path to x Feather file")
    parser.add_argument("--beta_dir", required=True, help="Beta files directory")
    parser.add_argument("--output", required=True, help="Output file base path")
    parser.add_argument("--use_single", action="store_true", help="Use float32")
    parser.add_argument("--ncores", type=int, default=4, help="Number of cores")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Configure GPU memory allocation
    configure_gpu_memory(args.ncores)
    
    # Start computation
    getf(args.x, args.beta_dir, args.use_single, args.output, args.ncores)

if __name__ == "__main__":
    # Set spawn method for cross-platform compatibility
    multiprocessing.set_start_method("spawn", force=True)
    print("Done")
    # Execute main function
    main()
