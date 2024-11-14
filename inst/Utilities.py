import tensorflow as tf

def is_gpu_available(print_status=True):
    """
    Checks if the system has any available GPUs for TensorFlow computations and 
    prints the number of GPUs available. It also returns a message indicating 
    whether the computation will be performed on GPU(s) or CPU.

    Returns:
        str: A message indicating if the system is using a GPU or CPU for computations.
    """
    # Get the list of available physical devices (GPUs) recognized by TensorFlow
    gpus = tf.config.list_physical_devices('GPU')
    
    # Print the number of GPUs available
    if print_status:
        num_gpus = len(gpus)
        print("Number of GPUs Available:", num_gpus)
        # Determine the compute device based on the presence of GPUs
        #if num_gpus > 0:
            #message = "The system is using GPU(s) for computations."
        #else:
            #message = "No GPU detected. The system will use the CPU for computations."
        #print(message)
    return len(gpus) > 0

def check_modules(module_list=None, print_status=True):
    """
    Checks the availability of specific Python modules.
    
    Args:
        module_list (list, optional): A list of module names to check for availability. If None, no action is taken.
        print_status (bool): If True, prints the availability of each module. Default is True.
    
    Returns:
        list: A list of module names that are not available (i.e., not installed).
        Returns None if no module list is provided.
    """
    # Check if module_list is provided; if not, print message and return None
    if not module_list:
        if print_status:
            print("No module list provided.")
        return None
    
    # List to store the names of unavailable modules
    unavailable_modules = []
    
    # Attempt to import each module and check if ImportError is raised
    for module in module_list:
        try:
            __import__(module)
        except ImportError:
            unavailable_modules.append(module)
    
    # Print the list of unavailable modules if print_status is True and any are missing
    if unavailable_modules and print_status:
        print("Modules not available:", ', '.join(unavailable_modules))
    
    return unavailable_modules
