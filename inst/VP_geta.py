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
#   `IASDT.R::VarPar_Compute`. The original functionality is extracted from the 
#   `Hmsc::computeVariancePartitioning` R function.
#   The implementation is optimized for parallel matrix multiplications using
#   TensorFlow. This script replaces the R function `geta`:
#
# geta = function(a){
#   switch(
#     class(hM$X)[1L],
#     matrix = {res = hM$X%*%(t(hM$Tr%*%t(a$Gamma)))},
#     list = {
#       res = matrix(NA,hM$ny,hM$ns)
#       for(j in 1:hM$ns) res[,j] =  hM$X[[j]]%*%(t(hM$Tr[j,]%*%t(a$Gamma)))
#     }
#   )
#   return(res)
# }
# 
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import multiprocessing
import argparse
import logging
import warnings
import os

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
import pyarrow as feather
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


def process_column(tr, x, gamma_chunk, col_offset, file_output):
    
    """
    Process a column chunk of the gamma matrix and compute x %*% (tr %*% gamma).

    Args:
        tr (tf.Tensor): Transformation matrix.
        x (tf.Tensor): Data matrix.
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
            
            tr_gamma = tf.reshape(gamma_chunk[:, col], [-1, 1])
            tr_gamma = tf.matmul(tr, tr_gamma, transpose_b=True)
            res = tf.matmul(x, tr_gamma, transpose_b=True)
            res = pd.DataFrame(res.numpy())
            
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
    
    for start in range(0, data.shape[1], chunk_size):
        yield data[:, start:start + chunk_size], start

# ======================================================================
# ======================================================================

def geta(file_tr, file_x, file_gamma, use_single, file_output, ncores, chunk_size):
    
    """
    Perform parallel matrix multiplication and save results.

    Args:
        file_tr (str): Path to the tr matrix Feather file.
        file_x (str): Path to the x matrix Feather file.
        file_Gamma (str): Path to the gamma matrix Feather file.
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
    
    tr = load_feather_tensor(file_tr, dtype)
    x = load_feather_tensor(file_x, dtype)
    gamma = load_feather_tensor(file_gamma, dtype)
    
    tasks = [(tr, x, chunk, offset, dtype, file_output) for chunk, offset in chunkify(gamma, chunk_size)]
    
    with get_reusable_executor(max_workers=ncores) as executor:
        futures = [executor.submit(process_column, *args) for args in tasks]
        
        for future in futures:
            try:
                future.result()
            except Exception as e:
                print(f"Job failed with exception: {e}")

# ======================================================================
# ======================================================================

def main():

    parser = argparse.ArgumentParser(
        description="Parallel TensorFlow-based matrix operations.")
    parser.add_argument(
        "--tr", required=True, help="Path to the tr Feather file")
    parser.add_argument(
        "--x", required=True, help="Path to the x Feather file")
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
    
    args = parser.parse_args()
    
    geta(args.tr, args.x, args.gamma, args.use_single, args.output, args.ncores, args.chunk_size)

if __name__ == "__main__":
    # Ensure this script can be run on both Windows and UNIX-based systems
    multiprocessing.set_start_method("spawn", force=True)
    print("Done")
    main()
