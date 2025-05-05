import os
import tensorflow as tf
import numpy as np
import logging
from concurrent.futures import ThreadPoolExecutor
from functools import partial

# ==============================================================================
# Configure logging
# ==============================================================================

# Set up logging to capture INFO-level messages for general runtime information
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ==============================================================================
# TensorFlow and Environment Configuration
# ==============================================================================

# Set TensorFlow to display only errors, suppressing warnings and info messages
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.get_logger().setLevel('ERROR')

# Disable oneDNN optimizations, which may add performance overhead for this use case
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

# ==============================================================================
# VarPar - Core Function for Matrix Computations with TensorFlow
# ==============================================================================

def VarPar(hM_X, hM_Tr, Post_Gamma, Post_Beta, use_single=False):
    
    """
    Perform a set of matrix multiplications with TensorFlow to calculate four matrices.

    Parameters:
    - hM_X (array-like): Input matrix X.
    - hM_Tr (array-like): Input matrix Tr.
    - Post_Gamma (array-like): Posterior matrix Gamma.
    - Post_Beta (array-like): Posterior matrix Beta.
    - use_single (bool): If True, uses tf.float32 (single precision); otherwise tf.float64.

    Returns:
    dict: Contains computed numpy arrays 'la', 'lf', 'lmu', 'lbeta'.
    """
    
    # Set the data type based on use_single flag (float32 or float64)
    dtype = tf.float32 if use_single else tf.float64

    # Convert input arrays to TensorFlow tensors of specified dtype
    hM_X = tf.convert_to_tensor(hM_X, dtype=dtype)
    hM_Tr = tf.convert_to_tensor(hM_Tr, dtype=dtype)
    Post_Gamma = tf.convert_to_tensor(Post_Gamma, dtype=dtype)
    Post_Beta = tf.convert_to_tensor(Post_Beta, dtype=dtype)

    # Ensure that hM_Tr and Post_Gamma are treated as 2D matrices
    if hM_Tr.shape[1] == 1:
        # Reshape to maintain 2D structure
        hM_Tr = tf.reshape(hM_Tr, [-1, 1])
    if Post_Gamma.shape[1] == 1:
        # Reshape to maintain 2D structure
        Post_Gamma = tf.reshape(Post_Gamma, [-1, 1])

    # Matrix operations
    la = tf.matmul(
        hM_X, tf.transpose(tf.matmul(hM_Tr, tf.transpose(Post_Gamma))))
    lf = tf.matmul(hM_X, Post_Beta)
    lmu = tf.transpose(tf.matmul(hM_Tr, tf.transpose(Post_Gamma)))
    lbeta = Post_Beta

    # Convert TensorFlow tensors to numpy arrays for R compatibility
    la = la.numpy()
    lf = lf.numpy()
    lmu = lmu.numpy()
    lbeta = lbeta.numpy()

    # Return the results as a dictionary of numpy arrays
    return {
        'la': la,
        'lf': lf,
        'lmu': lmu,
        'lbeta': lbeta
    }

# ==============================================================================
# process_post - Helper Function to Process Each Element in postList
# ==============================================================================

def process_post(post_data, hM_X, hM_Tr, use_single):
    
    """
    Helper function to call VarPar with data from a single element in postList.
    Handles any errors encountered during VarPar execution.

    Parameters:
    - post_data (dict): Dictionary containing 'Gamma' and 'Beta' matrices
    - hM_X (array-like): Matrix X for computations
    - hM_Tr (array-like): Matrix Tr for computations
    - use_single (bool): If True, uses float32; otherwise, uses float64

    Returns:
    dict or None: Results from VarPar, or None if an error occurred
    """
    
    try:
        # Call VarPar with extracted 'Gamma' and 'Beta' from post_data
        return VarPar(
            hM_X=hM_X,
            hM_Tr=hM_Tr,
            Post_Gamma=post_data['Gamma'],
            Post_Beta=post_data['Beta'],
            use_single=use_single)
    except Exception as e:
        # Log the exception and return None if an error occurs
        print("Exception in process_post:", e)
        return None

# ==============================================================================
# VarPar_Parallel - Execute VarPar in Parallel Across postList
# ==============================================================================

def VarPar_Parallel(hM_X, hM_Tr, postList, use_single=False, num_threads=None):
    
    """
    Executes VarPar for each element in postList in parallel using ThreadPoolExecutor.

    Parameters:
    - hM_X (array-like): Matrix X for computations
    - hM_Tr (array-like): Matrix Tr for computations
    - postList (list of dicts): List of dictionaries containing 'Gamma' and 'Beta'
    - use_single (bool): If True, uses float32; otherwise, uses float64
    - num_threads (int, optional): Number of threads to use. Defaults to CPU count.

    Returns:
    list: List of results from VarPar for each element in postList
    """
    
    # Set the number of threads based on num_threads parameter or system CPU count
    if num_threads is None:
        num_threads = os.cpu_count()
    
    # Create a partial function of process_post with preset arguments
    partial_process_post = partial(process_post, hM_X=hM_X, hM_Tr=hM_Tr, use_single=use_single)

    # Use ThreadPoolExecutor to run process_post in parallel across postList
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Map partial_process_post to each item in postList and collect results
        results = list(executor.map(partial_process_post, postList))

    return results

# ==============================================================================
# VarPar_Sequential - Execute VarPar Sequentially Across postList
# ==============================================================================

def VarPar_Sequential(hM_X, hM_Tr, postList, use_single=False):
    
    """
    Executes VarPar for each element in postList sequentially, not in parallel.

    Parameters:
    - hM_X (array-like): Matrix X for computations
    - hM_Tr (array-like): Matrix Tr for computations
    - postList (list of dicts): List of dictionaries containing 'Gamma' and 'Beta'
    - use_single (bool): If True, uses float32; otherwise, uses float64

    Returns:
    list: List of results from VarPar for each element in postList
    """

    # Iterate over each element in postList, calling process_post sequentially
    results = []
    for post_data in postList:
        result = process_post(post_data, hM_X, hM_Tr, use_single)
        results.append(result)
    return results

# ==============================================================================
# warmup
# ==============================================================================

def warmup():
    """A lightweight function to initialise necessary modules."""
    # Use os to list the current directory
    _ = os.listdir('.')
    
    # Create a small TensorFlow tensor
    tf.constant([0], dtype=tf.float32)
    
    # Perform a lightweight numpy operation
    _ = np.array([1, 2, 3]).mean()
    
    # Initialise logging
    logger = logging.getLogger("warmup_logger")
    logger.setLevel(logging.INFO)
    logger.info("Warm-up complete for logging")
    
    # Use ThreadPoolExecutor with a simple function
    def dummy_function(x):
        return x * 2

    with ThreadPoolExecutor(max_workers=2) as executor:
        _ = list(executor.map(partial(dummy_function), [1, 2]))
    
    # do not return objects
    pass
