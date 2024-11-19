import os
# Set TensorFlow logging level to show only errors
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# Disable oneDNN optimizations to prevent additional warnings
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

import logging
# disable Additional Logs Set the TensorFlow logger - show only critical errors
logging.getLogger('tensorflow').setLevel(logging.ERROR)

import tensorflow as tf
#tf.compat.v1.reset_default_graph()
#tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
# Reset all devices before making changes
tf.keras.backend.clear_session()
#tf.get_logger().setLevel(logging.ERROR)
#tf.get_logger().setLevel('ERROR')

import numpy as np
import rdata
import warnings # Suppress warnings during warm-up
# Suppress all TensorFlow warnings
warnings.filterwarnings("ignore", category=UserWarning, module="tensorflow")

import gc
import time
import argparse
#import json  # For formatting output as JSON if needed
#import sys
#print("Received arguments:", sys.argv)

# ======================================================================
# TensorFlow and Environment Configuration
# ======================================================================



# Use compatibility mode to avoid deprecated calls
#tf.compat.v1.reset_default_graph()  # This avoids the warning
#The name tf.reset_default_graph is deprecated. Please use tf.compat.v1.reset_default_graph instead."


# Flush the logger
#sys.stdout.flush()  # Ensure immediate output







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
# load_rds
# ======================================================================

def load_rds(file_path):
    
    """
    Load an RDS file using the rdata module.

    Parameters:
    - file_path (str): Path to the RDS file.

    Returns:
    - rds_data (dict): Loaded data from the RDS file.
    
    Raises:
    - FileNotFoundError: If the specified file does not exist.
    - ValueError: If there is an issue reading the RDS file or it is malformed.
    """

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The specified file does not exist: {file_path}")
    
    try:
        with open(file_path, 'rb') as f:
            rds_data = rdata.read_rds(f)
        return rds_data
    except ValueError as e:
        raise ValueError(f"Error reading RDS file at {file_path}: {e}")

# ======================================================================
# print_time
# ======================================================================

def print_time(message, verbose=True):
    
    """
    Print the current time along with a provided message.

    Parameters:
    - message (str): The message to print along with the current time.
    """
    
    if verbose:
        current_time = time.strftime("%I:%M:%S %p", time.localtime())
        print(f"{message} - {current_time}")

# ======================================================================
# crossprod_solve_old
# ======================================================================

# This is the old version of the function. This worked nicely on a small-scale
# models, but it needs to be optimized for larger data. See the new 
# version below.

def crossprod_solve_old(Dist1, Dist2, Denom, postEta, use_single=False):
    
    """
    Compute cross-product and solve matrix systems in TensorFlow.
    
    Parameters:
    - Dist1 (str or array-like): Distance matrix 1 or its file path.
    - Dist2 (str or array-like): Distance matrix 2 or its file path.
    - Denom (float): Denominator for exponent scaling.
    - postEta (str or array-like): 3D matrix (n x m x k) or its file path.
    - use_single (bool): If True, use float32 precision; otherwise, use float64.

    Returns:
    - results (list of np.ndarray): Flattened arrays containing results 
    for each submatrix in postEta.
    """
    
    # Set data type based on user input
    dtype = tf.float32 if use_single else tf.float64
    
    # Load and convert Dist1, Dist2, and postEta if they are file paths
    if isinstance(Dist1, str):
        Dist1 = load_rds(Dist1)
    if isinstance(Dist2, str):
        Dist2 = load_rds(Dist2)
    if isinstance(postEta, str):
        postEta = load_rds(postEta)
    
    # Check if loading was successful
    if Dist1 is None or Dist2 is None or postEta is None:
        print_time("Error: One or more input files could not be loaded.")
        return None
    
    # Convert distance matrices and denominator to TensorFlow tensors
    Dist1 = tf.convert_to_tensor(Dist1, dtype=dtype)
    Dist2 = tf.convert_to_tensor(Dist2, dtype=dtype)
    Denom_tensor = tf.convert_to_tensor(Denom, dtype=dtype)
        
    # Calculate the K matrices with element-wise exponential
    K11 = tf.exp(-Dist1 / Denom_tensor)
    K12 = tf.exp(-Dist2 / Denom_tensor)

    # Prepare an empty list to store results for each matrix in postEta
    results = []
    
    # Loop through each "slice" in the third dimension of postEta
    # (e.g., postEta[:, :, i])
    num_matrices = postEta.shape[2]  # Number of matrices (should be 32)
    
    for i in range(num_matrices):
        # Extract and convert the i-th 2D matrix from postEta
        Vect_tensor = tf.convert_to_tensor(postEta[:, :, i], dtype=dtype)

        # Attempt to solve the linear system for this slice
        try:
            # Solve K11 * x = Vect_tensor
            solved_vector = tf.linalg.solve(K11, Vect_tensor)
            # Compute cross-product with K12
            result = tf.matmul(K12, solved_vector, transpose_a=True)
            # Append the flattened result
            results.append(result.numpy().flatten())
        except tf.errors.InvalidArgumentError as e:
            print_time(f"Error solving for matrix {i} in postEta: {e}. Check dimensions and matrix properties.")
            # Append None if an error occurs for this slice
            results.append(None)
    
    # Return the list of results for all slices
    return results


# ======================================================================
# load_tensor_chunked
# ======================================================================

def load_tensor_chunked(file_or_array, dtype, chunk_size=1000, threshold_mb=2000):
    
    """
    Convert file or array to a TensorFlow tensor in chunks to handle large data.
    
    Parameters:
    - file_or_array (str or array-like): Input data or file path.
    - dtype (tf.DType): TensorFlow data type (e.g., tf.float64).
    - chunk_size (int): Size of chunks to process at a time.
    - threshold_mb (int): Threshold for chunking in megabytes (default 2000 MB).

    Returns:
    - tensor (tf.Tensor): Combined tensor from chunks.
    """

    
    if isinstance(file_or_array, str):
        file_or_array = rdata.read_rds(file_or_array)

    # Ensure positive integer chunk size
    chunk_size = max(1, int(chunk_size))
    threshold_bytes = threshold_mb * 1024**2

    if isinstance(file_or_array, str):
        file_or_array = load_rds(file_or_array)
    
    # If the array is too large, split into chunks
    # 2 GB threshold
    if isinstance(file_or_array, np.ndarray) and file_or_array.nbytes > threshold_bytes:
        chunks = []
        for i in range(0, file_or_array.shape[0], chunk_size):
            chunk = file_or_array[i : i + chunk_size]
            chunks.append(tf.convert_to_tensor(chunk, dtype=dtype))
        return tf.concat(chunks, axis=0)
    else:
        return tf.convert_to_tensor(file_or_array, dtype=dtype)

# ======================================================================
# compute_k_matrices
# ======================================================================

@tf.function
def compute_k_matrices(Dist1, Dist2, Denom_tensor):
    
    """
    Compute K matrices as element-wise exponentials.
    
    Parameters:
    - Dist1 (tf.Tensor): Distance matrix 1.
    - Dist2 (tf.Tensor): Distance matrix 2.
    - Denom_tensor (tf.Tensor): Denominator for scaling.

    Returns:
    - K11, K12 (tf.Tensor, tf.Tensor): Computed matrices.
    """
    
    return tf.exp(-Dist1 / Denom_tensor), tf.exp(-Dist2 / Denom_tensor)

# ======================================================================
# crossprod_solve_old2 - faster than crossprod_solve2, but sill slow
# ======================================================================

def crossprod_solve_old2(Dist1, Dist2, Denom, postEta, use_single=False, save=False, path_out=None, chunk_size=1000, threshold_mb=2000):
    
    """
    Compute cross-product and solve matrix systems in TensorFlow.
    Optionally save the output to an RDS file.

    Parameters:
    - Dist1 (str or array-like): Distance matrix 1 or its file path.
    - Dist2 (str or array-like): Distance matrix 2 or its file path.
    - Denom (float): Denominator for exponent scaling.
    - postEta (str or array-like): 3D matrix (n x m x k) or its file path.
    - use_single (bool): If True, use float32 precision; otherwise, use float64.
    - save (bool): If True, save results to an RDS file.
    - path_out (str): File path to save the results if `save` is True.
    - chunk_size (int): Size of chunks for loading large matrices.

    Returns:
    - results (list of np.ndarray): Flattened arrays containing results for 
        each submatrix in postEta.
    """

    # Ensure chunk_size is an integer
    chunk_size = int(chunk_size)

    # Set data type based on user input
    dtype = tf.float32 if use_single else tf.float64

    # Load and convert input matrices
    Dist1 = load_tensor_chunked(Dist1, dtype, chunk_size=chunk_size, threshold_mb = threshold_mb)
    Dist2 = load_tensor_chunked(Dist2, dtype, chunk_size=chunk_size, threshold_mb = threshold_mb)
    if isinstance(postEta, str):
        postEta = load_tensor_chunked(load_rds(postEta), dtype, chunk_size=chunk_size, threshold_mb = threshold_mb)

    # Check if loading was successful
    if Dist1 is None or Dist2 is None or postEta is None:
        print_time("Error: One or more input files could not be loaded.")
        return None

    # Convert Denom to a TensorFlow tensor
    Denom_tensor = tf.convert_to_tensor(Denom, dtype=dtype)

    # Calculate the K matrices
    K11, K12 = compute_k_matrices(Dist1, Dist2, Denom_tensor)
    
    # Free memory
    del Dist1, Dist2, Denom_tensor
    gc.collect()

    # Reshape postEta for batch processing
    num_matrices = postEta.shape[2]
    postEta_reshaped = tf.reshape(postEta, [postEta.shape[0], postEta.shape[1] * num_matrices])

    # Free memory
    del postEta
    gc.collect()

    # Solve linear systems in batch
    results = tf.linalg.solve(K11, postEta_reshaped)

    # Compute cross-products in batch
    results = tf.matmul(K12, results, transpose_a=True)

    # Reshape results back to the original structure
    results = tf.reshape(results, [num_matrices, -1])

    # Convert results to numpy arrays
    results = [res.numpy().flatten() for res in tf.split(results, num_matrices)]

    # Save to RDS if requested
    if save:
        if not path_out:
            raise ValueError("A valid path_out must be provided when save=True.")
        try:
            with open(path_out, 'wb') as f:
                rdata.write_rds(results, f)
        except Exception as e:
            print_time(f"Error saving results to {path_out}: {e}")

    return results

# ======================================================================
# compute_pairwise_distances
# ======================================================================

def compute_pairwise_distances(s1_chunk, s2, dtype=tf.float64):
    
    """
    Compute pairwise distances between rows of s1_chunk and s2.

    Parameters:
    - s1_chunk (tensor): A 2D tensor of coordinates (n_chunk x d).
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
    Compute pairwise distances between rows of s1 and s2 in chunks to manage memory usage.

    Parameters:
    - s1 (tensor): A 2D tensor of coordinates (n x d).
    - s2 (tensor): Optional; A 2D tensor of coordinates (m x d) to compute distances against.
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
# crossprod_solve
# ======================================================================

def crossprod_solve(s1, s2, denom, postEta, use_single=False, save=False, path_out=None, chunk_size=1000, threshold_mb=2000, verbose=False):
    
    """
    Compute cross-product and solve matrix systems in TensorFlow.
    Optionally save the output to an RDS file.

    Parameters:
    - s1 (str or array-like): Distance matrix 1 or its file path.
    - s2 (str or array-like): Distance matrix 2 or its file path.
    - denom (float): Denominator for exponent scaling.
    - postEta (str or array-like): 3D matrix (n x m x k) or its file path.
    - use_single (bool): If True, use float32 precision; otherwise, use float64.
    - save (bool): If True, save results to an RDS file.
    - path_out (str): File path to save the results if save is True.
    - chunk_size (int): Size of chunks for loading large matrices.
    - verbose (bool): If True, print execution time messages; otherwise, suppress them.
    
    Returns:
    - results (list of np.ndarray): Flattened arrays containing results for each submatrix in postEta.
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
        s1_tensor = load_rds(s1)
        print_time("  >>  Convert to tensor", verbose)
        s1_tensor = tf.convert_to_tensor(s1_tensor, dtype=dtype)
        print_time("  >>  Compute d11", verbose)
        Dist1 = compute_distances_chunked(
            s1=s1_tensor, dtype=dtype, chunk_size=chunk_size)

    print_time("Processing coordinate matrix - s2", verbose)
    if isinstance(s2, str):
        print_time("  >>  Reading data", verbose)
        s2_tensor = load_rds(s2)
        print_time("  >>  Convert to tensor", verbose)
        s2_tensor = tf.convert_to_tensor(s2_tensor, dtype=dtype)
        print_time("  >>  Compute d12", verbose)
        Dist2 = compute_distances_chunked(
            s1=s1_tensor, s2=s2_tensor, dtype=dtype, chunk_size=chunk_size)
    
    print_time("  >>  Free memory", verbose)
    del s1_tensor, s2_tensor, s1, s2
    gc.collect()
    
    # |||||||||||||||||||||||||||||||||||
    # Load and convert postEta
    # |||||||||||||||||||||||||||||||||||
    print_time("Load and convert postEta", verbose)
    #if isinstance(postEta, str):
        #postEta = load_tensor_chunked(postEta, dtype, chunk_size, threshold_mb)

    if isinstance(postEta, str):
        try:
            postEta = load_tensor_chunked(postEta, dtype, chunk_size, threshold_mb)
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
    print_time("  >>  k12", verbose)
    K12 = tf.exp(-Dist2 / Denom_tensor)
    
    # Free memory
    print_time("Free memory", verbose)
    del Dist1, Dist2, Denom_tensor, denom
    gc.collect()

    # |||||||||||||||||||||||||||||||||||
    # Reshape postEta for batch processing
    # |||||||||||||||||||||||||||||||||||
    print_time("Reshape postEta for batch processing", verbose)
    num_matrices = postEta.shape[2]
    postEta_reshaped = tf.reshape(postEta, [postEta.shape[0], postEta.shape[1] * num_matrices])
    
    # Free memory
    print_time("Free memory", verbose)
    del postEta
    gc.collect()

    # |||||||||||||||||||||||||||||||||||
    # Solve linear systems in batch
    # |||||||||||||||||||||||||||||||||||
    print_time("Solve linear systems in batch", verbose)
    results = tf.linalg.solve(K11, postEta_reshaped)

    # |||||||||||||||||||||||||||||||||||
    # Compute cross-products in batch
    # |||||||||||||||||||||||||||||||||||
    print_time("Compute cross-products in batch", verbose)
    results = tf.matmul(K12, results, transpose_a=True)

    # |||||||||||||||||||||||||||||||||||
    # Reshape results back to the original structure
    # |||||||||||||||||||||||||||||||||||
    print_time("Reshape results back to the original structure", verbose)
    results = tf.reshape(results, [num_matrices, -1])

    # Free memory
    print_time("Free memory", verbose)
    del postEta_reshaped, K11, K12
    gc.collect()

    # |||||||||||||||||||||||||||||||||||
    # Convert results to numpy arrays
    # |||||||||||||||||||||||||||||||||||
    print_time("Convert results to numpy arrays", verbose)
    results = [res.numpy().flatten() for res in tf.split(results, num_matrices)]
    results = [arr.flatten().tolist() for arr in results]

    # |||||||||||||||||||||||||||||||||||
    # Save to RDS if requested
    # |||||||||||||||||||||||||||||||||||

    if save:
        if not path_out:
            raise ValueError("A valid path_out must be provided when save=True.")
        try:
            # Flatten arrays before saving
            results2 = [arr[0] for arr in results]
            rdata.write_rds(path_out, results2)
        except Exception as e:
            print_time(f"Error saving results to {path_out}: {e}", verbose)

    # |||||||||||||||||||||||||||||||||||
    # Print processing time
    # |||||||||||||||||||||||||||||||||||
    elapsed_time = time.time() - start_time
    print_time(f"Finished in {elapsed_time:.1f} seconds", verbose)

    return results

# ==============================================================================
# warmup
# ==============================================================================

def warmup():
    """A lightweight function to initialize necessary modules and resources."""
    start_time = time.time()  # Record the start time
    
    # Use os to get the current working directory
    _ = os.getcwd()
    
    # Create a small TensorFlow tensor
    tf.constant([0], dtype=tf.float32)
    
    # Perform a lightweight numpy operation
    _ = np.array([0, 1, 2]).sum()
    
    # Suppress rdata warnings during warm-up
    dummy_rds = b'x\x9cK\xccK.\xca\xcf+I\xd5Q(\xc9H-.QHI,IT\xc8/.-\xceHL\xb1\x00\xb2\xd1\x19\xb4'
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            _ = rdata.parser.parse_data(dummy_rds)
        except Exception:
            # Catch exceptions since we're not using a real RDS file
            pass
    
    # Perform garbage collection
    gc.collect()
    
    # Log the elapsed time
    elapsed_time = time.time() - start_time
    del elapsed_time
    
    # Do not return objects
    pass

# ==============================================================================
# ==============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run crossprod_solve function.")
    parser.add_argument("--s1", type=str, required=True, help="Path to s1 file.")
    parser.add_argument("--s2", type=str, required=True, help="Path to s2 file.")
    parser.add_argument("--postEta", type=str, required=True, help="Path to postEta file.")
    parser.add_argument("--path_out", type=str, required=True, help="Path to output file.")
    parser.add_argument("--denom", type=float, required=True, help="Denominator value.")
    parser.add_argument("--chunk_size", type=int, default=1000, help="Chunk size.")
    parser.add_argument("--threshold_mb", type=int, default=2000, help="Threshold in MB.")
    parser.add_argument("--use_single", action="store_true", help="Use single precision.")
    parser.add_argument("--save", action="store_true", help="Save the output.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output.")
    args = parser.parse_args()

    try:
        result = crossprod_solve(
            s1=args.s1,
            s2=args.s2,
            denom=args.denom,
            postEta=args.postEta,
            use_single=args.use_single,
            save=args.save,
            path_out=args.path_out,
            chunk_size=args.chunk_size,
            threshold_mb=args.threshold_mb,
            verbose=args.verbose,
        )

        if result is None:
            raise ValueError("crossprod_solve returned None.")
        print("Done")
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
