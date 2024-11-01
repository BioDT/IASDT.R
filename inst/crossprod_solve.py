import os
import tensorflow as tf
import numpy as np
import rdata
import xarray as xr
import pandas as pd


# =======================================================
# TensorFlow and Environment Configuration
# =======================================================

# Set TensorFlow logging level to show only errors
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.get_logger().setLevel('ERROR')

# Disable oneDNN optimizations to prevent additional warnings
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

# check GPU availability
# print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))


# =======================================================
# load_rds
# =======================================================

def load_rds(file_path):
    """
    Load an RDS file using the rdata module.

    Parameters:
    - file_path (str): Path to the RDS file.

    Returns:
    - rds_data (dict): Loaded data from the RDS file.
    
    Raises:
    - FileNotFoundError: If the specified file does not exist.
    - rdata.errors.ReadError: If there is an issue reading the RDS file.
    """
    
    try:
        with open(file_path, 'rb') as f:
            rds_data = rdata.read_rds(f)
        return rds_data
    except (FileNotFoundError, rdata.errors.ReadError) as e:
        print(f"Error loading RDS file at {file_path}: {e}")
        return None


# =======================================================
# crossprod_solve
# =======================================================

def crossprod_solve(Dist1, Dist2, Denom, List, use_single=False):
    
    """
    Compute cross-product and solve matrix systems in TensorFlow.
    
    Parameters:
    - Dist1 (str or array-like): Distance matrix 1 or its file path.
    - Dist2 (str or array-like): Distance matrix 2 or its file path.
    - Denom (float): Denominator for exponent scaling.
    - List (str or array-like): 3D matrix (n x m x k) or its file path.
    - use_single (bool): If True, use float32 precision; otherwise, use float64.

    Returns:
    - results (list of np.ndarray): Flattened arrays containing results for each submatrix in List.
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
    
    # Loop through each "slice" in the third dimension of List (e.g., List[:, :, i])
    num_matrices = List.shape[2]  # Number of matrices (should be 32)
    
    
    for i in range(num_matrices):
        # Extract and convert the i-th 2D matrix from List
        Vect_tensor = tf.convert_to_tensor(List[:, :, i], dtype=dtype)

        # Attempt to solve the linear system for this slice
        try:
            # Solve K11 * x = Vect_tensor
            solved_vector = tf.linalg.solve(K11, Vect_tensor)
            # Compute cross-product with K12
            result = tf.matmul(K12, solved_vector, transpose_a=True)
            # Append the flattened result
            results.append(result.numpy().flatten())
        except tf.errors.InvalidArgumentError as e:
            print(f"Error solving for matrix {i} in List: {e}. Check dimensions and matrix properties.")
            # Append None if an error occurs for this slice
            results.append(None)
    
    # Return the list of results for all slices
    return results
