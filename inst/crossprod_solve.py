import os
# Set TensorFlow logging level to show only errors
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# Disable oneDNN optimizations to prevent additional warnings
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

import logging
# Disable Additional Logs Set the TensorFlow logger - show only critical errors
# https://stackoverflow.com/questions/55081911/
logging.getLogger('tensorflow').setLevel(logging.ERROR)

import tensorflow as tf
# Reset all devices before making changes
tf.keras.backend.clear_session()
# show only errors, no warnings
# https://www.tensorflow.org/api_docs/python/tf/get_logger
tf.get_logger().setLevel('ERROR')

# https://www.tensorflow.org/api_docs/python/tf/compat/v1/reset_default_graph
# tf.compat.v1.reset_default_graph()
# tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

import warnings
# Suppress all TensorFlow warnings
warnings.filterwarnings("ignore", category=UserWarning, module="tensorflow")

import numpy as np
import pandas as pd
import gc
import sys
import time
import argparse
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

def load_tensor_chunked(file_or_array, dtype, chunk_size=1000, threshold_mb=2000):
    
    """
    Convert file or array to a TensorFlow tensor in chunks for large data handling.

    Parameters:
    - file_or_array (str or array-like): Input data or file path.
    - dtype (tf.DType): TensorFlow data type (e.g., tf.float64).
    - chunk_size (int, optional): Size of chunks to process at a time.
    - threshold_mb (int, optional): Threshold for chunking in megabytes (default 2000 MB).

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

def compute_pairwise_distances(s1_chunk, s2, dtype=tf.float64):
    
    """
    Compute Euclidean distances between rows of s1_chunk and s2.

    Parameters:
    - s1_chunk: A 2D tensor of shape (n, 2), where n is the number of rows.
    - s2 (tensor): A 2D tensor of coordinates (m x d).
    - dtype (tf.DType): Data type for the computation.

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

def compute_distances_chunked(s1, s2=None, dtype=tf.float64, chunk_size=1000):
    
    """
    Compute Euclidean distances between rows of s1 and s2 in chunks to manage memory usage.

    Parameters:
    - s1 (tensor): A 2D tensor of coordinates (n x d).
    - s2 (tensor): Optional; A 2D tensor of coordinates (m x d) to compute distances against. If not provided, `s1` is used.
    - dtype (tf.DType): Data type for computations.
    - chunk_size (int): Number of rows to process in each chunk to avoid memory overload.

    Returns:
    - Dist (tensor): A 2D tensor containing the computed pairwise distances.
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
        chunk_distances = compute_pairwise_distances(s1_chunk, s2, dtype)
        dist_list.append(chunk_distances)

    # Concatenate all chunks to form the full distance matrix
    Dist = tf.concat(dist_list, axis=0)
    return Dist

# ======================================================================
# solve_and_multiply
# ======================================================================

def solve_and_multiply(K11, K12, postEta):
    """
    Solve K11 * X = postEta for each column of postEta, extracting the first 
    item of each list, then compute tf.matmul(K12, X, transpose_a=True) for 
    each column.

    Parameters:
    - K11: TensorFlow tensor of shape [n, n].
    - K12: TensorFlow tensor of shape [n, p].
    - postEta: TensorFlow tensor of shape [n, m].

    Returns:
    - TensorFlow tensor: A 2D tensor of shape (p, m) for combined result after tf.matmul for all columns of postEta.
    """
    
    num_columns = postEta.shape[1]  # Number of columns in postEta
    final_results = []  # To store the results of tf.matmul for each column

    for col in range(num_columns):
        # Extract the column as a 1D tensor (shape: [n])
        postEta_col = postEta[:, col]

        # Expand dimensions to match input requirements for tf.linalg.solve
        # Shape: [n, 1]
        postEta_col_expanded = tf.expand_dims(postEta_col, axis=1)

        # Solve K11 * X = postEta_col
        # Shape: [n, 1]
        X = tf.linalg.solve(K11, postEta_col_expanded)

        # Perform tf.matmul(K12, X, transpose_a=True)
        # Shape: [p, 1]
        result = tf.matmul(K12, X, transpose_a=True)

        # Append the result
        final_results.append(result)

        # Clean up memory
        del postEta_col, postEta_col_expanded, X, result
        gc.collect()

    # Merge all results into a single tensor
    # Shape: [p, m]
    final_results = tf.concat(final_results, axis=1)

    return final_results

# ======================================================================
# crossprod_solve
# ======================================================================

def crossprod_solve(s1, s2, denom, postEta, use_single=False, path_out=None, chunk_size=1000, threshold_mb=2000, verbose=False):
    
    """
    Compute cross-product and solve matrix systems in TensorFlow and save the 
    output to an feather file.

    Parameters:
    - s1 (str or array-like): Distance matrix 1 or its file path.
    - s2 (str or array-like): Distance matrix 2 or its file path.
    - denom (float): Denominator for exponent scaling.
    - postEta (str or array-like): 3D matrix (n x m x k) or its file path.
    - use_single (bool): If True, use float32 precision; otherwise, use float64.
    - path_out (str): File path to save the results.
    - chunk_size (int): Size of chunks for loading large matrices.
    - threshold_mb: Threshold in megabytes for chunking.
    - verbose (bool): If True, print execution time messages; otherwise, suppress them.
    
    Returns:
    - path_out: File path to save the results.
    """
    
    print_time(f"Starting processing (TensorFlow v{tf.__version__})", verbose)

    if gpus:
        print_time(f"  >>  GPUs detected: {len(gpus)}", verbose)
    else:
        print_time("  >>  No GPUs detected; using CPU.", verbose)

    # Record the start time
    start_time = time.time()
        
    # Ensure chunk_size is an integer
    chunk_size = int(chunk_size)

    # Set data type based on user input
    dtype = tf.float32 if use_single else tf.float64

    # |||||||||||||||||||||||||||||||||||
    # Reading coordinate matrices 
    # |||||||||||||||||||||||||||||||||||
    print_time("Processing coordinate matrix - s1", verbose)
    if isinstance(s1, str):
        print_time("  >>  Reading data", verbose)
        s1_tensor = load_feather(s1)
        print_time("  >>  Convert to tensor", verbose)
        s1_tensor = tf.convert_to_tensor(s1_tensor, dtype=dtype)
        gc.collect()
        print_time("  >>  Compute d11", verbose)
        Dist1 = compute_distances_chunked(s1=s1_tensor, dtype=dtype, chunk_size=chunk_size)

    print_time("Processing coordinate matrix - s2", verbose)
    if isinstance(s2, str):
        print_time("  >>  Reading data", verbose)
        s2_tensor = load_feather(s2)
        print_time("  >>  Convert to tensor", verbose)
        s2_tensor = tf.convert_to_tensor(s2_tensor, dtype=dtype)
        print_time("  >>  Compute d12", verbose)
        Dist2 = compute_distances_chunked(s1=s1_tensor, s2=s2_tensor, dtype=dtype, chunk_size=chunk_size)
    
    print_time("  >>  Free memory", verbose)
    del s1_tensor, s2_tensor, s1, s2
    gc.collect()
    
    # |||||||||||||||||||||||||||||||||||
    # Load and convert postEta
    # |||||||||||||||||||||||||||||||||||
    print_time("Load and convert postEta", verbose)
    if isinstance(postEta, str):
        try:
            postEta = load_feather(postEta)
            postEta  = load_tensor_chunked(postEta , dtype, chunk_size, threshold_mb)
        except Exception as e:
            print_time(f"Error loading postEta: {e}", verbose)
            sys.exit(1)
    
    # |||||||||||||||||||||||||||||||||||
    # Convert denom to a TensorFlow tensor
    # |||||||||||||||||||||||||||||||||||
    print_time("Convert denom to a TensorFlow tensor", verbose)
    Denom_tensor = tf.convert_to_tensor(denom, dtype=dtype)
    
    # |||||||||||||||||||||||||||||||||||
    # Check if data was loaded successfully
    # |||||||||||||||||||||||||||||||||||
    print_time("Check if data was loaded successfully", verbose)
    if Dist1 is None or Dist2 is None or postEta is None or Denom_tensor is None:
        print_time("Error: One or more input files could not be loaded.", verbose)
        return None

    # |||||||||||||||||||||||||||||||||||
    # Calculate the K matrices
    # |||||||||||||||||||||||||||||||||||
    print_time("Calculate the K matrices", verbose)
    
    print_time("  >>  k11", verbose)
    K11 = tf.exp(-Dist1 / Denom_tensor)
    del Dist1
    
    print_time("  >>  k12", verbose)
    K12 = tf.exp(-Dist2 / Denom_tensor)
    # Retrieve the shape as a list [rows, columns]
    K12_shape = K12.shape.as_list()
    # Dynamically ensure the tensor's shape
    K12 = tf.ensure_shape(K12, K12_shape)

    # Free memory
    print_time("  >>  Free memory", verbose)
    del Dist2, Denom_tensor, denom
    gc.collect()

    # |||||||||||||||||||||||||||||||||||
    # Solve & cross-product
    # |||||||||||||||||||||||||||||||||||
    print_time("Solve & cross-product", verbose)
    results = solve_and_multiply(K11, K12, postEta)
    
    # Free memory
    print_time("  >>  Free memory", verbose)
    del K11, K12, postEta
    gc.collect()
    
    # |||||||||||||||||||||||||||||||||||
    # Convert results to numpy arrays
    # |||||||||||||||||||||||||||||||||||
    print_time("Convert results to pandas dataframe", verbose)
    results = pd.DataFrame(results)
    
    # |||||||||||||||||||||||||||||||||||
    # Save to feather
    # |||||||||||||||||||||||||||||||||||
    print_time("Save to feather", verbose)

    if not path_out:
        raise ValueError("A valid path_out must be provided.")
    try:
        results.to_feather(path_out)
    except Exception as e:
        print_time(f"Error saving results to {path_out}: {e}", verbose)

    # |||||||||||||||||||||||||||||||||||
    # Print processing time
    # |||||||||||||||||||||||||||||||||||
    elapsed_time = time.time() - start_time
    print_time(f"Finished in {elapsed_time:.1f} seconds", verbose)

    return path_out

# ==============================================================================
# ==============================================================================

# Allow `crossprod_solve` function from the command line (system / system2 in R).

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run crossprod_solve function.")

    parser = argparse.ArgumentParser(
    description=(
        "Run crossprod_solve to compute cross-product and solve linear systems "
        "using TensorFlow, and save the output in Feather format."))
    parser.add_argument("--s1", type=str, required=True, help="Path to s1 file.")
    parser.add_argument("--s2", type=str, required=True, help="Path to s2 file.")
    parser.add_argument("--postEta", type=str, required=True, help="Path to postEta file.")
    parser.add_argument("--path_out", type=str, required=True, help="Path to output file.")
    parser.add_argument("--denom", type=float, required=True, help="Denominator value.")
    parser.add_argument("--chunk_size", type=int, default=1000, help="Chunk size for processing data in memory.")
    parser.add_argument("--threshold_mb", type=int, default=2000, help="Memory threshold in MB for chunking large data.")
    parser.add_argument("--use_single", action="store_true", help="Use single precision.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output.")
    args = parser.parse_args()

    try:
        result = crossprod_solve(
            s1=args.s1,
            s2=args.s2,
            denom=args.denom,
            postEta=args.postEta,
            use_single=args.use_single,
            path_out=args.path_out,
            chunk_size=args.chunk_size,
            threshold_mb=args.threshold_mb,
            verbose=args.verbose)

        if result is None:
            raise ValueError("crossprod_solve returned None.")
        print("Done")
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
