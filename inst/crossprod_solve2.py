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

# =======================================================

def load_rds(file_path):
    """
    Load an RDS file using the rdata module.
    """
    
    try:
        with open(file_path, 'rb') as f:
            rds_data = rdata.read_rds(f)
        return rds_data
    except (FileNotFoundError, rdata.errors.ReadError) as e:
        print(f"Error loading RDS file at {file_path}: {e}")
        return None

# =======================================================

@tf.function
def load_tensor(file_or_array, dtype):
    """
    Convert file or array to a TensorFlow tensor.
    """
    if isinstance(file_or_array, str):
        file_or_array = load_rds(file_or_array)
    return tf.convert_to_tensor(file_or_array, dtype=dtype)


# =======================================================

@tf.function
def compute_k_matrices(Dist1, Dist2, Denom_tensor):
    """
    Compute K matrices as element-wise exponentials.
    """
    return tf.exp(-Dist1 / Denom_tensor), tf.exp(-Dist2 / Denom_tensor)

# =======================================================

def crossprod_solve(Dist1, Dist2, Denom, List, use_single=False):
    """
    Compute cross-product and solve matrix systems in TensorFlow.
    """
    # Set data type based on user input
    dtype = tf.float32 if use_single else tf.float64

    # Load and convert input matrices
    Dist1 = load_tensor(Dist1, dtype)
    Dist2 = load_tensor(Dist2, dtype)
    if isinstance(List, str):
        List = load_rds(List)

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

    # Solve linear systems in batch
    solved_vectors = tf.linalg.solve(K11, List_reshaped)

    # Compute cross-products in batch
    results = tf.matmul(K12, solved_vectors, transpose_a=True)

    # Reshape results back to the original structure
    results = tf.reshape(results, [num_matrices, -1])

    # Convert results to numpy arrays for further processing if needed
    return [res.numpy().flatten() for res in tf.split(results, num_matrices)]
