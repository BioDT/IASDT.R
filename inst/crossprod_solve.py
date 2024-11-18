import os
import tensorflow as tf
import numpy as np
import rdata
import warnings # Suppress warnings during warm-up
import gc

# ======================================================================
# TensorFlow and Environment Configuration
# ======================================================================

# Set TensorFlow logging level to show only errors
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.get_logger().setLevel('ERROR')

# Disable oneDNN optimizations to prevent additional warnings
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

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
# crossprod_solve_old
# ======================================================================

# This is the old version of the function. This worked nicely on a small-scale
# models, but it needs to be optimized for larger data. See the new 
# version below.

def crossprod_solve_old(Dist1, Dist2, Denom, List, use_single=False):
    
    """
    Compute cross-product and solve matrix systems in TensorFlow.
    
    Parameters:
    - Dist1 (str or array-like): Distance matrix 1 or its file path.
    - Dist2 (str or array-like): Distance matrix 2 or its file path.
    - Denom (float): Denominator for exponent scaling.
    - List (str or array-like): 3D matrix (n x m x k) or its file path.
    - use_single (bool): If True, use float32 precision; otherwise, use float64.

    Returns:
    - results (list of np.ndarray): Flattened arrays containing results 
    for each submatrix in List.
    """
    
    # Set data type based on user input
    dtype = tf.float32 if use_single else tf.float64
    
    # Load and convert Dist1, Dist2, and List if they are file paths
    if isinstance(Dist1, str):
        Dist1 = load_rds(Dist1)
    if isinstance(Dist2, str):
        Dist2 = load_rds(Dist2)
    if isinstance(List, str):
        List = load_rds(List)
    
    # Check if loading was successful
    if Dist1 is None or Dist2 is None or List is None:
        print("Error: One or more input files could not be loaded.")
        return None
    
    # Convert distance matrices and denominator to TensorFlow tensors
    Dist1 = tf.convert_to_tensor(Dist1, dtype=dtype)
    Dist2 = tf.convert_to_tensor(Dist2, dtype=dtype)
    Denom_tensor = tf.convert_to_tensor(Denom, dtype=dtype)
        
    # Calculate the K matrices with element-wise exponential
    K11 = tf.exp(-Dist1 / Denom_tensor)
    K12 = tf.exp(-Dist2 / Denom_tensor)

    # Prepare an empty list to store results for each matrix in List
    results = []
    
    # Loop through each "slice" in the third dimension of List
    # (e.g., List[:, :, i])
    num_matrices = List.shape[2]  # Number of matrices (should be 32)
    
    for i in range(num_matrices):
        # Extract and convert the i-th 2D matrix from List
        Vect_tensor = tf.convert_to_tensor(List[:, :, i], dtype=dtype)

        # Attempt to solve the linear system for this slice
        try:
            # Solve K11 * x = Vect_tensor
            result = tf.linalg.solve(K11, Vect_tensor)
            # Compute cross-product with K12
            result = tf.matmul(K12, result, transpose_a=True)
            # Append the flattened result
            results.append(result.numpy().flatten())
        except tf.errors.InvalidArgumentError as e:
            print(f"Error solving for matrix {i} in List: {e}. Check dimensions and matrix properties.")
            # Append None if an error occurs for this slice
            results.append(None)
    
    # Return the list of results for all slices
    return results

# ======================================================================
# load_tensor
# ======================================================================

@tf.function
def load_tensor(file_or_array, dtype):
    """
    Convert file or array to a TensorFlow tensor.
    """
    if isinstance(file_or_array, str):
        file_or_array = load_rds(file_or_array)
    tensor = tf.convert_to_tensor(file_or_array, dtype=dtype)
    return tensor

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
    
    threshold_bytes = threshold_mb * 1024**2
    
    if isinstance(file_or_array, str):
        file_or_array = rdata.read_rds(file_or_array)

    #Ensure chunk_size is an integer
    chunk_size = int(chunk_size)

    # If the array is too large, split into chunks
    if isinstance(file_or_array, np.ndarray) and file_or_array.nbytes > threshold_bytes
        chunks = []
        for i in range(0, file_or_array.shape[0], chunk_size):
            chunk = file_or_array[i : i + chunk_size]
            chunks.append(tf.convert_to_tensor(chunk, dtype=dtype))
            del chunk    
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
# crossprod_solve
# ======================================================================

def crossprod_solve(Dist1, Dist2, Denom, List, use_single=False, save=False, file_path=None, chunk_size=1000, threshold_mb=2000):
    """
    Compute cross-product and solve matrix systems in TensorFlow.
    Optionally save the output to an RDS file.

    Parameters:
    - Dist1 (str or array-like): Distance matrix 1 or its file path.
    - Dist2 (str or array-like): Distance matrix 2 or its file path.
    - Denom (float): Denominator for exponent scaling.
    - List (str or array-like): 3D matrix (n x m x k) or its file path.
    - use_single (bool): If True, use float32 precision; otherwise, use float64.
    - save (bool): If True, save results to an RDS file.
    - file_path (str): File path to save the results if `save` is True.
    - chunk_size (int): Size of chunks for loading large matrices.

    Returns:
    - results (list of np.ndarray): Flattened arrays containing results for each submatrix in List.
    """

    # Ensure chunk_size is an integer
    chunk_size = int(chunk_size)

    # Set data type based on user input
    dtype = tf.float32 if use_single else tf.float64

    # Load and convert input matrices
    Dist1 = load_tensor_chunked(Dist1, dtype, chunk_size=chunk_size, threshold_mb = threshold_mb)
    Dist2 = load_tensor_chunked(Dist2, dtype, chunk_size=chunk_size, threshold_mb = threshold_mb)
    if isinstance(List, str):
        List = load_tensor_chunked(load_rds(List), dtype, chunk_size=chunk_size, threshold_mb = threshold_mb)

    # Check if loading was successful
    if Dist1 is None or Dist2 is None or List is None:
        print("Error: One or more input files could not be loaded.")
        return None

    # Convert Denom to a TensorFlow tensor
    Denom_tensor = tf.convert_to_tensor(Denom, dtype=dtype)

    # Calculate the K matrices
    K11, K12 = compute_k_matrices(Dist1, Dist2, Denom_tensor)

    # Reshape List for batch processing
    num_matrices = List.shape[2]
    List_reshaped = tf.reshape(List, [List.shape[0], List.shape[1] * num_matrices])

    # Free memory
    del Dist1, Dist2, Denom_tensor List
    gc.collect()

    # Solve linear systems in batch
    solved_vectors = tf.linalg.solve(K11, List_reshaped)
    
    del K11, List_reshaped
    gc.collect()

    # Compute cross-products in batch
    results = tf.matmul(K12, solved_vectors, transpose_a=True)

    # Reshape results back to the original structure
    results = tf.reshape(results, [num_matrices, -1])

    # Convert results to numpy arrays
    results = [res.numpy().flatten() for res in tf.split(results, num_matrices)]

    # Save to RDS if requested
    if save:
        if not file_path:
            raise ValueError("A valid file_path must be provided when save=True.")
        try:
            with open(file_path, 'wb') as f:
                rdata.write_rds(results, f)
        except Exception as e:
            print(f"Error saving results to {file_path}: {e}")

    return results

# ==============================================================================
# warmup
# ==============================================================================

def warmup():
    """A lightweight function to initialize necessary modules."""
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
    
    # do not return objects
    pass
