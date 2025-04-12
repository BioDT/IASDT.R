# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# - Authors: Ahmed El-Gabbas [ahmed.el-gabbas@ufz.de]
# - Last updated: 27.11.2024
# - Affiliation: Helmholtz Centre for Environmental Research - UFZ, Germany
# - License: MIT License
#
# - Description:
#   
#   This Python script aids in computing variance partitioning for Hmsc models.
#   It is designed to be called within the R function 
#   `IASDT.R::variance_partitioning_compute`. The original functionality is 
#   extracted from the  `Hmsc::computeVariancePartitioning` R function.
#   The implementation is optimized for parallel matrix multiplications using 
#   TensorFlow. This script replaces the R function `gemu`:
#
# gemu = function(a){
#   res = t(hM$Tr%*%t(a$Gamma))
#   return(res)
# }
# 
#   CITATION: El-Gabbas, A. (2025). IASDT.R: Modelling the distribution of invasive alien plant species in Europe. https://doi.org/10.5281/zenodo.1483438, R package version 0.1.X; https://biodt.github.io/IASDT.R/.
#
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import multiprocessing
import argparse
import logging
import warnings
import os
import traceback

# Disable Additional Logs Set the TensorFlow logger - show only critical errors
# https://stackoverflow.com/questions/55081911/
logging.getLogger('tensorflow').setLevel(logging.ERROR)

# Set TensorFlow logging level to show only errors
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# Disable oneDNN optimizations to prevent additional warnings
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

# Suppress all TensorFlow warnings
warnings.filterwarnings("ignore", category=UserWarning, module="tensorflow")

import tensorflow as tf
# Reset all devices before making changes
tf.keras.backend.clear_session()
# show only errors, no warnings
# https://www.tensorflow.org/api_docs/python/tf/get_logger
tf.get_logger().setLevel('ERROR')

# https://www.tensorflow.org/api_docs/python/tf/compat/v1/reset_default_graph
# tf.compat.v1.reset_default_graph()
# tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

import numpy as np
import pandas as pd
import pyarrow.feather as feather
from loky import get_reusable_executor

# ======================================================================
# TensorFlow and Environment Configuration
# ======================================================================

gpus = tf.config.experimental.list_physical_devices('GPU')

if gpus:
    try:
        # Enable memory growth
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
# ======================================================================

def load_feather_tensor(file_path, dtype):
    """
    Load a Feather file and return its data as a TensorFlow tensor.

    Args:
        file_path (str): Path to the Feather file.
        dtype (tf.DType): TensorFlow data type (e.g., tf.float32).

    Returns:
        tf.Tensor: Data loaded as a TensorFlow tensor.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If Feather data reading fails.
    """
        
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: The specified file '{file_path}' does not exist.")
    try:
        feather_data = feather.read_feather(file_path).to_numpy()
        # Ensure data is contiguous in memory
        feather_data = np.ascontiguousarray(feather_data)
        feather_data = tf.convert_to_tensor(feather_data, dtype=dtype)
        return feather_data
    except ValueError as e:
        raise ValueError(f"Failed to read Feather file at {file_path}: {e}")

# ======================================================================
# ======================================================================

def process_column(tr, gamma_chunk, col_offset, file_output):
    
    """
    Process a column chunk of the gamma matrix.

    Args:
        tr (tf.Tensor): Transformation matrix.
        gamma_chunk (tf.Tensor): Chunk of the gamma matrix.
        col_offset (int): Column offset for indexing.
        file_output (str): Output file path.
    """

    try:
        for col in range(gamma_chunk.shape[1]):
            
            col_index = col + col_offset
            
            # Dynamically generate file names - create a 4-digit padded index
            file_index = f"{col_index + 1:04d}"
            feather_file = f"{file_output.replace('.feather', f'_{file_index}.feather')}"

            if os.path.exists(feather_file):
                continue
                        
            # Ensure it's a column vector
            gamma = tf.reshape(gamma_chunk[:, col], [-1, 1])
            
            # Matrix multiplication: tr %*% t(gamma)
            tr_gamma = tf.matmul(tr, gamma, transpose_b=True)
            
            # convert to pandas DataFrame
            res = pd.DataFrame(tf.transpose(tr_gamma).numpy())
            
            # Save as Feather
            res.to_feather(feather_file)
        
    except Exception as e:
        print(f"Job failed with exception: {e}")
        raise

# ======================================================================
# ======================================================================

def chunkify(data, chunk_size):
    """
    Split data into chunks of specified size.

    Args:
        data (tf.Tensor): Input data to be chunked.
        chunk_size (int): Size of each chunk.

    Yields:
        tuple: Chunk of data and its start index.
    """

    # Iterate over the data tensor
    for start in range(0, data.shape[1], chunk_size):
        # Yield the chunk and its start index
        yield data[:, start:start + chunk_size], start

# ======================================================================
# ======================================================================

def gemu(file_tr, file_gamma, use_single, file_output, ncores, chunk_size):
    
    """
    Perform parallel matrix multiplication and save results.

    Args:
        file_tr (str): Path to the tr matrix Feather file.
        file_gamma (str): Path to the gamma matrix Feather file.
        use_single (bool): Whether to use single precision (float32).
        file_output (str): Path to the output file.
        ncores (int): Number of CPU cores for parallel processing.
        chunk_size (int): Size of chunks for processing.
    """

    if chunk_size <= 0:
        raise ValueError("Chunk size must be a positive integer.")
    
    if ncores <= 0:
        raise ValueError("Number of cores must be a positive integer.")

    dtype = tf.float32 if use_single else tf.float64
    
    # Load the tr and gamma matrices
    tr = load_feather_tensor(file_tr, dtype)
    gamma = load_feather_tensor(file_gamma, dtype)
    
    # Split the gamma matrix into chunks
    tasks = [(tr, chunk, offset, file_output) for chunk, offset in chunkify(gamma, chunk_size)]
    
    # Process each chunk in parallel
    with get_reusable_executor(max_workers=ncores, timeout=600) as executor:
        futures = []
        for args in tasks:
            try:
                futures.append(executor.submit(process_column, *args))
            except Exception as e:
                print(f"Error submitting task: {e}")

        for future in futures:
            try:
                future.result()
            except Exception as e:
                print(f"Job failed with exception: {e}")
                traceback.print_exc()

# ======================================================================
# ======================================================================

def main():

    # Create an argument parser
    parser = argparse.ArgumentParser(description="Parallel TensorFlow-based matrix operations.")
    parser.add_argument(
        "--tr", required=True, help="Path to the tr Feather file")
    parser.add_argument(
        "--gamma", required=True, help="Path to the gamma Feather file")
    parser.add_argument(
        "--output", required=True, help="Output Feather file path")
    parser.add_argument(
        "--use_single", action="store_true", 
        help="Use single precision for computations.")
    parser.add_argument(
        "--ncores", type=int, default=4, help="Number of CPU cores")
    parser.add_argument(
        "--chunk_size", type=int, default=50, 
        help="Size of chunks for processing")

    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Call the gemu function with the provided arguments
    gemu(args.tr, args.gamma, args.use_single, args.output, args.ncores, args.chunk_size)

# ======================================================================
# ======================================================================

if __name__ == "__main__":
    # Ensure this script can be run on both Windows and UNIX-based systems
    multiprocessing.set_start_method("spawn", force=True)
    print("Done")
    main()
