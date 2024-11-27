# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# - Authors: Ahmed El-Gabbas [ahmed.el-gabbas@ufz.de]
# - Last updated: 27.11.2024
# - Affiliation: Helmholtz Centre for Environmental Research - UFZ, Germany
# - License: MIT License
#
# - Description:
#
#   This Python script aids in computing predictions for the latent factor of 
#   Hmsc models. It is designed to be called within the R function 
#   `IASDT.R::predictLF`. The original functionality is extracted from the 
#   `Hmsc::predictLatentFactor` R function.
#   The implementation is optimized for parallel matrix multiplications using #   TensorFlow. This script was tested for single spatial random effects using 
#   the GPP method.
# 
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import os
import logging
import warnings
import time
import argparse
import contextlib
import gc
import sys

# Set TensorFlow logging level to show only errors
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# Disable oneDNN optimizations to prevent additional warnings
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

# Disable Additional Logs Set the TensorFlow logger - show only critical errors
# https://stackoverflow.com/questions/55081911/
logging.getLogger('tensorflow').setLevel(logging.ERROR)

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
# file_path
# ======================================================================

def load_feather(file_path):
    
    """
    Load a Feather file and return its data as a NumPy array.

    Parameters:
    - file_path (str): Path to the Feather file.

    Returns:
    - np.ndarray: Data loaded from the Feather file.

    Raises:
    - FileNotFoundError: If the file does not exist.
    - ValueError: If the Feather file cannot be read.
    """

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found:: {file_path}")
    try:
        feather_data = feather.read_feather(file_path).to_numpy()
        return feather_data
    except ValueError as e:
        raise ValueError(f"Failed to read Feather file at {file_path}: {e}")

# ======================================================================
# print_time
# ======================================================================

def print_time(message, verbose=True):
    
    """
    Print the current time along with a provided message  (optional).

    Parameters:
    - message (str): The message to print along with the current time.
    - verbose (bool, optional): If True, print the message. Defaults to True.
    """

    if verbose:
        current_time = time.strftime("%I:%M:%S %p", time.localtime())
        print(f"{message} - {current_time}")

# ======================================================================
# load_tensor_chunked
# ======================================================================

def load_tensor_chunked(
        file_or_array, dtype, chunk_size=1000, threshold_mb=2000):
    
    """
    Convert file or array to a TensorFlow tensor in chunks for large 
    data handling.

    Parameters:
    - file_or_array (str or array-like): Input data or file path.
    - dtype (tf.DType): TensorFlow data type (e.g., tf.float64).
    - chunk_size (int, optional): Size of chunks to process at a time.
    - threshold_mb (int, optional): Threshold for chunking in megabytes 
    (default 2000 MB).

    Returns:
    - tensor (tf.Tensor): Combined tensor from chunks.
    """

    if isinstance(file_or_array, str):
        #file_or_array = rdata.read_rds(file_or_array)
        file_or_array = load_feather(file_or_array)
    
    # Ensure positive integer chunk size
    chunk_size = max(1, int(chunk_size))
    threshold_bytes = threshold_mb * 1024**2

    # If the array is too large, split into chunks - 2 GB threshold
    if isinstance(file_or_array, np.ndarray) and file_or_array.nbytes > threshold_bytes:
        chunks = []
        for i in range(0, file_or_array.shape[0], chunk_size):
            chunk = file_or_array[i : i + chunk_size]
            chunks.append(tf.convert_to_tensor(chunk, dtype=dtype))
        return tf.concat(chunks, axis=0)
    else:
        return tf.convert_to_tensor(file_or_array, dtype=dtype)

# ======================================================================
# compute_pairwise_distances
# ======================================================================

def compute_pairwise_distances(s1_chunk, s2):
    
    """
    Compute Euclidean distances between rows of s1_chunk and s2.

    Parameters:
    - s1_chunk: A 2D tensor of shape (n, 2), where n is the number of rows.
    - s2 (tensor): A 2D tensor of coordinates (m x d).

    Returns:
    - dist_matrix (tensor): A 2D tensor of pairwise distances (n_chunk x m).
    """

    # Compute squared Euclidean distance
    squared_diff = tf.expand_dims(s1_chunk, 1) - tf.expand_dims(s2, 0)
    squared_dist = tf.reduce_sum(tf.square(squared_diff), axis=-1)
    
    # Take square root to get the Euclidean distance
    dist_matrix = tf.sqrt(squared_dist)

    return dist_matrix

# ======================================================================
# compute_distances_chunked
# ======================================================================

def compute_distances_chunked(
        s1, s2=None, dtype=tf.float64, chunk_size=1000):
    
    """
    Compute Euclidean distances between rows of s1 and s2 in chunks to 
    manage memory usage.

    Parameters:
    - s1 (tensor): A 2D tensor of coordinates (n x d).
    - s2 (tensor): Optional; A 2D tensor of coordinates (m x d) to compute 
    distances against. If not provided, `s1` is used.
    - dtype (tf.DType): Data type for computations.
    - chunk_size (int): Number of rows to process in each chunk to avoid 
    memory overload.

    Returns:
    - dist (tensor): A 2D tensor containing the computed pairwise distances.
    """
    
    # If s2 is not provided, use s1 itself (for self-distances)
    if s2 is None:
        s2 = s1

    # Ensure s1 and s2 are tensors
    s1 = tf.convert_to_tensor(s1, dtype=dtype)
    s2 = tf.convert_to_tensor(s2, dtype=dtype)
    
    # Initialize an empty list to store results
    dist_list = []

    # Iterate over chunks of s1
    for start_idx in range(0, s1.shape[0], chunk_size):
        end_idx = min(start_idx + chunk_size, s1.shape[0])
        s1_chunk = s1[start_idx:end_idx]

        # For each chunk of s1, compute the distances to all of s2
        chunk_distances = compute_pairwise_distances(s1_chunk, s2)
        dist_list.append(chunk_distances)

    # Concatenate all chunks to form the full distance matrix
    dist = tf.concat(dist_list, axis=0)
    return dist

# ======================================================================
# solve_and_multiply
# ======================================================================

# The function supports chunked processing (solve_chunk_size > 1) for solving 
# the linear system and cross-product computations, reducing memory usage for 
# large datasets.
# If chunking is not used (solve_chunk_size == 1), the function processes each 
# column of post_eta sequentially.
# Memory Optimization: Explicitly cleans up memory with gc.collect() after 
# processing chunks to reduce memory overhead.
# Dynamic Logging: Logs progress dynamically to track the exact stages 
# of computation.

def solve_and_multiply(k11, k12, post_eta, log_fn=None, solve_chunk_size=50):
    
    """
    Solve the linear system and compute cross-products for each column of 
    post_eta efficiently with optional chunking.
    
    Solve k11 * x = post_eta
    Parameters:
    - k11 (Tensor): Kernel matrix of shape [n, n].
    - k12 (Tensor): Kernel matrix of shape [n, p].
    - post_eta (Tensor): Matrix of shape [n, m] to solve for each column.
    - log_fn (function): Logging function for progress messages.
    - solve_chunk_size (int): Number of columns to process in chunks. 
    Use 1 for single-column processing.
    
    If solve_chunk_size == 1 (no chunking), yhe function processes one column 
    of post_eta at a time, ensuring isolated computation for each column.
    If solve_chunk_size > 1 (chunking), columns of post_eta are processed in 
    groups (chunks), where each chunk contains multiple columns. This reduces 
    the number of calls to tf.linalg.solve and tf.matmul, improving efficiency.
    
    Returns:
    - Tensor: Results of the cross-product computations, shape [p, m].
    """
    
    # Number of columns in post_eta
    num_columns = post_eta.shape[1]
    # List to store results for each chunk or column
    final_results = []
    
    # Chunked processing
    if solve_chunk_size > 1:
        if log_fn:
            log_fn("\n  >>  Processing solve_and_multiply in chunks")

        for start_col in range(0, num_columns, solve_chunk_size):
            end_col = min(start_col + solve_chunk_size, num_columns)

            # Extract chunk of columns from post_eta
            post_eta_chunk = post_eta[:, start_col:end_col]
            if log_fn:
                log_fn(f"  >>  >>  Processing columns {start_col + 1} to {end_col} of {num_columns}")

            # Solve k11 * x = post_eta_chunk
            x_chunk = tf.linalg.solve(k11, post_eta_chunk)

            # Compute cross-product: tf.matmul(k12, x_chunk, transpose_a=True)
            chunk_result = tf.matmul(k12, x_chunk, transpose_a=True)

            # Append chunk results
            final_results.append(chunk_result)

            # Clean up memory
            del post_eta_chunk, x_chunk, chunk_result
            gc.collect()

        # Combine all chunks into a single tensor
        final_results = tf.concat(final_results, axis=1)


    # Single-column processing
    else:
        if log_fn:
            log_fn("\n  >>  Processing solve_and_multiply without chunking\n")

        for col in range(num_columns):
            if log_fn:
                log_fn(f"  >>  Processing column {col + 1} of {num_columns}")

            # Extract a single column from post_eta
            post_eta_col = tf.expand_dims(post_eta[:, col], axis=1)

            # Solve k11 * x = post_eta_col
            x = tf.linalg.solve(k11, post_eta_col)

            # Compute cross-product: tf.matmul(k12, x, transpose_a=True)
            results = tf.matmul(k12, x, transpose_a=True)

            # Append results
            final_results.append(results)

            # Clean up memory
            del post_eta_col, x, results
            gc.collect()

        # Combine all columns into a single tensor
        final_results = tf.concat(final_results, axis=1)

    return final_results

# ======================================================================
# crossprod_solve
# ======================================================================

# The crossprod_solve function performs the following tasks:
# - Read and Process Large Matrices: Reads data from Feather files for two 
# distance matrices (s1, s2) and a target matrix (post_eta).
# - Computes pairwise Euclidean distances in chunks for efficient memory 
# management.
# - Matrix Operations: Constructs kernel matrices (k11, k12) using the 
# computed distances and a scaling factor (denom).
# - Solves the linear system K11⋅X=postEtaK11⋅X=post_eta and computes 
# cross-products K12T⋅XK12T⋅X, handling chunking for large datasets.
# Output:
# - Saves the results as a Feather file.
# - Logs the progress and runtime in a .log file.
        
def crossprod_solve(
        s1, s2, denom, post_eta, use_single=False, path_out=None, 
        chunk_size=1000, threshold_mb=2000, verbose=True, solve_chunk_size=50):
    
    """
    Compute cross-product and solve matrix systems in TensorFlow, saving 
    # results to a Feather file.

    Workflow:
    1. Load data (`s1`, `s2`, `post_eta`) from Feather files.
    2. Compute pairwise distances to create kernel matrices (`k11`, `k12`).
    3. Solve the equation `k11 * x = post_eta` for each column (or chunk) of post_eta.
    4. Compute cross-products `k12^T * x` for each column (or chunk) of post_eta.
    5. Save results to a Feather file.

    Parameters:
    - s1 (str or array-like): Path to or data for the first distance matrix.
    - s2 (str or array-like): Path to or data for the second distance matrix.
    - denom (float): Denominator for scaling the exponent in kernel matrices.
    - post_eta (str or array-like): Path to or data for the target matrix.
    - use_single (bool): Use single precision (`float32`) if True, else double (`float64`).
    - path_out (str): File path to save results.
    - chunk_size (int): Chunk size for pairwise distance computations.
    - threshold_mb (int): Memory threshold in MB for chunking.
    - verbose (bool): If True, print progress to the log file.
    - solve_chunk_size (int): Number of columns to process simultaneously in solve_and_multiply.

    Returns:
    - str: Path to the saved Feather file.
    """
    
    
    log_file = path_out.replace(".feather", ".log")
    
    with open(log_file, "a") as log:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            
            # Print to file with time stamp
            def log_and_flush(msg, verbose=True):
                if verbose:
                    print_time(msg)
                    log.flush()
            
            # Print to file without time stamp
            def log_and_flush2(msg, verbose=True):
                if verbose:
                    print(msg)
                    log.flush()
            
            # Record the start time
            start_time = time.time()
            
            # Log environment info
            log_and_flush2("\n", verbose)
            log_and_flush2("=" * 80, verbose)
            log_and_flush("Starting crossprod_solve", verbose)
            log_and_flush2("=" * 80, verbose)
            
            # Log environment info
            log_and_flush2("\nEnvironment Info:", verbose)
            log_and_flush2(f"    Current Working Directory: {os.getcwd()}", verbose)
            log_and_flush2(f"    Python Version: {sys.version}", verbose)
            log_and_flush2(f"    TensorFlow Version: {tf.__version__}", verbose)
            log_and_flush2(f"    OS: {os.name}", verbose)
            log_and_flush2(f"    Platform: {sys.platform}", verbose)

            # Log user inputs
            log_and_flush2("\nUser Inputs:", verbose)
            log_and_flush2(f"    s1: {s1}", verbose)
            log_and_flush2(f"    s2: {s2}", verbose)
            log_and_flush2(f"    denom: {denom}", verbose)
            log_and_flush2(f"    post_eta: {post_eta}", verbose)
            log_and_flush2(f"    use_single: {use_single}", verbose)
            log_and_flush2(f"    path_out: {path_out}", verbose)
            log_and_flush2(f"    chunk_size: {chunk_size}", verbose)
            log_and_flush2(f"    threshold_mb: {threshold_mb}", verbose)
            log_and_flush2(f"    solve_chunk_size: {solve_chunk_size}", verbose)
            
            # Check for GPUs
            log_and_flush2("\n\nChecking GPU", verbose)
            if gpus:
                log_and_flush2(f"  >>  GPUs detected: {len(gpus)}\n", verbose)
            else:
                log_and_flush2("  >>  No GPUs detected; using CPU.\n", verbose)
        
            # Ensure chunk_size is an integer
            chunk_size = int(chunk_size)
        
            # Set data type
            dtype = tf.float32 if use_single else tf.float64
        
            # |||||||||||||||||||||||||||||||||||
            # Processing s1
            # |||||||||||||||||||||||||||||||||||
            log_and_flush("Processing s1", verbose)
                
            try:
                s1_tensor = load_feather(s1) if isinstance(s1, str) else s1
                s1_tensor = tf.convert_to_tensor(s1_tensor, dtype=dtype)
                gc.collect()
                dist1 = compute_distances_chunked(
                    s1_tensor, dtype=dtype, chunk_size=chunk_size)
            except Exception as e:
                log_and_flush(f"Error processing s1: {e}", verbose)
                raise
        
            # |||||||||||||||||||||||||||||||||||
            # Processing s2
            # |||||||||||||||||||||||||||||||||||
        
            log_and_flush("Processing s2", verbose)
            
            try:
                s2_tensor = load_feather(s2) if isinstance(s2, str) else s2
                s2_tensor = tf.convert_to_tensor(s2_tensor, dtype=dtype)
                dist2 = compute_distances_chunked(
                    s1_tensor, s2_tensor, dtype=dtype, chunk_size=chunk_size)
            except Exception as e:
                log_and_flush(f"Error processing s2: {e}", verbose)
                raise
        
            del s1_tensor, s2_tensor
            gc.collect()
            
            # |||||||||||||||||||||||||||||||||||
            # Processing post_eta
            # |||||||||||||||||||||||||||||||||||
            
            log_and_flush("Processing post_eta", verbose)
            
            try:
                post_eta = load_feather(post_eta) if isinstance(post_eta, str) else post_eta
                post_eta = load_tensor_chunked(
                    post_eta, dtype=dtype, chunk_size=chunk_size, 
                    threshold_mb=threshold_mb)
            except Exception as e:
                log_and_flush(f"Error loading post_eta: {e}", verbose)
                raise
            
            # |||||||||||||||||||||||||||||||||||
            # Compute k11 and k12
            # |||||||||||||||||||||||||||||||||||
            
            log_and_flush("Calculating K matrices", verbose)
                        
            try:
                # Convert denom to a TensorFlow tensor
                denom_tensor = tf.convert_to_tensor(denom, dtype=dtype)
                
                log_and_flush("Calculating k11", verbose)
                k11 = tf.exp(-dist1 / denom_tensor)
                del dist1

                log_and_flush("Calculating k12", verbose)
                k12 = tf.exp(-dist2 / denom_tensor)
                # Retrieve the shape as a list [rows, columns]
                #K12_shape = k12.shape.as_list()
                # Dynamically ensure the tensor's shape
                #k12 = tf.ensure_shape(k12, K12_shape)
                del dist2, denom_tensor
                gc.collect()

            except Exception as e:
                log_and_flush(f"Error calculating K matrices: {e}", verbose)
                raise
            
            # |||||||||||||||||||||||||||||||||||
            # Solve and multiply
            # |||||||||||||||||||||||||||||||||||
            log_and_flush("Solving and computing cross-products", verbose)
            
            try:
                results = solve_and_multiply(
                    k11, k12, post_eta, log_fn=log_and_flush,
                    solve_chunk_size=solve_chunk_size)
                del k11, k12, post_eta
                gc.collect()
            except Exception as e:
                log_and_flush(f"Error in solve_and_multiply: {e}", verbose)
                raise

            # |||||||||||||||||||||||||||||||||||
            # Save results
            # |||||||||||||||||||||||||||||||||||
            
            log_and_flush("Saving results to Feather file", verbose)
            
            try:
                pd.DataFrame(results.numpy()).to_feather(path_out)
            except Exception as e:
                log_and_flush(f"Error saving results: {e}", verbose)
                raise
            
            # |||||||||||||||||||||||||||||||||||
            # Print processing time
            # |||||||||||||||||||||||||||||||||||
            
            elapsed_time = time.time() - start_time
            
            log_and_flush(f"Elapsed time: {time.strftime('%H:%M:%S', time.gmtime(elapsed_time))}", verbose)
            log_and_flush2("=" * 80, verbose)
            log_and_flush2("\n\n", verbose)

    return path_out

# ==============================================================================
# ==============================================================================

# Allow `crossprod_solve` function from the command line (system / system2 in R).

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description=(
    "Run crossprod_solve to compute cross-product and solve linear systems "
    "using TensorFlow, and save the output in Feather format."))
    
    parser.add_argument(
        "--s1", type=str, required=True,
        help="Path to the first distance matrix (s1).")
    parser.add_argument(
        "--s2", type=str, required=True, 
        help="Path to the second distance matrix (s2).")
    parser.add_argument(
        "--post_eta", type=str, required=True, 
        help="Path to the target matrix (post_eta).")
    parser.add_argument(
        "--path_out", type=str, required=True, 
        help="Path to save the output Feather file.")
    parser.add_argument(
        "--denom", type=float, required=True, 
        help="Scaling denominator for kernel matrices.")
    parser.add_argument(
        "--chunk_size", type=int, default=1000, 
        help="Chunk size for distance computations.")
    parser.add_argument(
        "--threshold_mb", type=int, default=2000, 
        help="Memory threshold in MB for chunking.")
    parser.add_argument(
        "--use_single", action="store_true", help="Use single precision for computations.")
    parser.add_argument(
        "--verbose", action="store_true", help="Enable verbose output.")
    parser.add_argument(
        "--solve_chunk_size", type=int, default=50, 
        help="Chunk size for solve_and_multiply.")

    args = parser.parse_args()
    
    try:
        result = crossprod_solve(
            s1=args.s1,
            s2=args.s2,
            denom=args.denom,
            post_eta=args.post_eta,
            use_single=args.use_single,
            path_out=args.path_out,
            chunk_size=args.chunk_size,
            threshold_mb=args.threshold_mb,
            verbose=args.verbose,
            solve_chunk_size=args.solve_chunk_size)

        if result is None:
            raise ValueError("crossprod_solve returned None.")
        print("Done")
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
