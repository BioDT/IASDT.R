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
#   The implementation is optimized for parallel matrix multiplications using #   TensorFlow. This script replaces the R function `getf`:
#
# getf = function(a){
#   switch(
#     class(hM$X)[1L],
#     matrix = {res = hM$X%*%(a$Beta)},
#     list = {
#       res = matrix(NA,hM$ny,hM$ns)
#       for(j in 1:hM$ns) res[,j] = hM$X[[j]]%*%a$Beta[,j]
#     }
#   )
#   return(res)
# }
# 
#   CITATION: El-Gabbas, A. (2025). IASDT.R: Modelling the distribution of invasive alien plant species in Europe. 10.5281/zenodo.14954739, R package version 0.1.X; https://biodt.github.io/IASDT.R/.
# 
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import multiprocessing
import argparse
import logging
import warnings
import os
import re
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

# Set a fixed memory limit (e.g., 64 GB) to avoid exceeding the GPU memory

def configure_gpu_memory(ncores, total_memory_mb=64 * 1024):
    """
    Configure GPU memory limits dynamically based on the number of cores.
    
    Args:
        ncores (int): Number of cores for parallel processing.
        total_memory_mb (int): Total available GPU memory in MB.
    """
    gpus = tf.config.list_physical_devices('GPU')
    
    if gpus:
        try:
            for gpu in gpus:
                # Validate ncores
                if ncores > 0:
                    memory_limit_per_worker = total_memory_mb // ncores
                else:
                    raise ValueError("Number of cores (ncores) must be greater than zero.")
                
                # Apply memory limit
                tf.config.set_logical_device_configuration(
                    gpu,
                    [tf.config.LogicalDeviceConfiguration(memory_limit=memory_limit_per_worker)]
                )
                # Enable memory growth
                #tf.config.experimental.set_memory_growth(gpu, True)
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

def process_beta_file(x, beta_file, task_index, dtype, file_output):
    
    """
    Perform matrix multiplication for a Beta chunk and save the results.

    Args:
        x (tf.Tensor): Input matrix x.
        beta_file (str): Path to the Beta chunk Feather file.
        dtype (tf.DType): Data type for calculations.
        file_output (str): Path to the output Feather file.

    Raises:
        ValueError: If matrix dimensions are incompatible.
    """
        
    try:
        # Dynamically generate file names
        file_index = f"{task_index + 1:04d}"
        feather_file = file_output.replace(".feather", f"_{file_index}.feather")
        
        # Skip processing if the output file already exists
        if os.path.exists(feather_file):
            return
        
        # Load the beta chunk
        beta_chunk = load_feather_tensor(beta_file, dtype)

        # Check dimensions for compatibility
        if beta_chunk.shape[0] != x.shape[1]:
            raise ValueError(f"Matrix size-incompatible: x columns ({x.shape[1]}) must match Beta rows ({beta_chunk.shape[0]})")

        # Retry logic for failed processing
        max_retries = 5
        for attempt in range(1, max_retries + 1):
            try:
                # Perform matrix multiplication: x %*% Beta
                res = tf.matmul(x, beta_chunk)
                
                # Save result as Feather
                pd.DataFrame(res.numpy()).to_feather(feather_file)
                # Exit the loop if successful
                return

            except Exception as e:
                if attempt == max_retries:
                    print(f"Max retries reached for {feather_file}. Skipping.")
                    raise

    except Exception as e:
        print(f"Job failed with exception: {e}")
        raise

# ======================================================================
# ======================================================================

def validate_and_get_beta_files(beta_dir, prefix="VP_Beta_", extension=".feather"):

    """
    Validate and retrieve files matching the specified prefix and extension in a directory.
    
    This function ensures that all files in the directory matching the pattern 
    '{prefix}{suffix}{extension}' are present, with a continuous range of suffix indices. 
    If files are missing, it raises an error.

    Args:
        beta_dir (str): Path to the directory containing Beta files.
        prefix (str): Prefix for file matching (default is "VP_Beta_").
        extension (str): File extension for matching (default is ".feather").
    
    Returns:
        list: Sorted list of paths to the validated Beta files.

    Raises:
        FileNotFoundError: If no matching files are found or if there are missing files in the range.
    """

    # Compile the pattern to match files (e.g., "VP_Beta_0001.feather")
    pattern = re.compile(rf"^{re.escape(prefix)}(\d+){re.escape(extension)}$")
    beta_files = []
    file_indices = []

    # Iterate through the files in the directory
    for f in os.listdir(beta_dir):
        match = pattern.match(f)
        if match:
            beta_files.append(os.path.join(beta_dir, f))
            # Extract numeric suffix as an integer
            file_indices.append(int(match.group(1)))

    # Check if any files matched
    if not beta_files:
        raise FileNotFoundError(f"No files matching the pattern '{prefix}*{extension}' found in {beta_dir}")
    
    # Validate the range of file indices
    file_indices = sorted(file_indices)
    # Full expected range
    expected_indices = list(range(file_indices[0], file_indices[-1] + 1))
    missing_files = set(expected_indices) - set(file_indices)

    if missing_files:
        missing_files_str = ", ".join([f"{prefix}{i:04d}{extension}" for i in sorted(missing_files)])
        raise FileNotFoundError(f"The following files are missing: {missing_files_str}")

    # Sort the list of beta files to ensure consistent order
    return sorted(beta_files)

# ======================================================================
# ======================================================================

def getf(file_x, beta_dir, use_single, file_output, ncores):
    
    """
    Perform parallel matrix multiplication and save results.

    Args:
        file_x (str): Path to the x matrix Feather file.
        beta_dir (str): Directory containing the Beta Feather files.
        use_single (bool): Whether to use single precision (float32).
        file_output (str): Path to the output file.
        ncores (int): Number of CPU cores for parallel processing.
    """
    
    dtype = tf.float32 if use_single else tf.float64
    x = load_feather_tensor(file_x, dtype)
    
    # Get list of Beta files
    beta_files = validate_and_get_beta_files(beta_dir)

    try:
        if ncores == 1:
            for task_index, beta_file in enumerate(beta_files):
                process_beta_file(x, beta_file, task_index, dtype, file_output)
        else:
            with get_reusable_executor(max_workers=ncores, timeout=600) as executor:
                futures = [
                    executor.submit(process_beta_file, x, beta_file, task_index, dtype, file_output)
                    for task_index, beta_file in enumerate(beta_files)
                ]

                for future in futures:
                    try:
                        future.result()
                    except Exception as e:
                        print(f"Worker failed with exception: {e}")
                        traceback.print_exc()
        
    except Exception as e:
        print(f"Error in parallel processing: {e}")
        raise

# ======================================================================
# ======================================================================

def main():

    parser = argparse.ArgumentParser(
        description="Parallel TensorFlow-based matrix operations.")
    parser.add_argument(
        "--x", required=True, help="Path to the x Feather file")
    parser.add_argument(
        "--beta_dir", required=True, 
        help="Directory containing the Beta Feather files")
    parser.add_argument(
        "--output", required=True, help="Output Feather file path")
    parser.add_argument(
        "--use_single", action="store_true", 
        help="Use single precision for computations.")
    parser.add_argument(
        "--ncores", type=int, default=4, help="Number of CPU cores")

    args = parser.parse_args()
    
    # Configure GPU memory
    configure_gpu_memory(args.ncores)

    getf(args.x, args.beta_dir, args.use_single, args.output, args.ncores)

if __name__ == "__main__":
    # Ensure this script can be run on both Windows and UNIX-based systems
    multiprocessing.set_start_method("spawn", force=True)
    print("Done")
    main()
