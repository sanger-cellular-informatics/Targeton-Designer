import logging
import sys, os


def logs_to_log(module_name: str) -> None:

    """
    Function logs_to_log() will log out all of current logging output to a log file,
    which sits under "/src/logs" directory. It will grab the logs outputted from a file
    in which it's instantiated is added.

    This function take one argument module_name and returns None.  
    For example: "ParentModule.Module". 

    Here, default logger level is DEBUG. 
    This means you can log INFO, ERROR, WARNING and EXCEPTION.
    """


    # Get directory for logs
    dir_path = os.path.dirname(os.path.realpath(__file__))
    parent_dir = os.path.dirname(dir_path)
    log_dir = os.path.join(parent_dir, 'logs/')
    module_name = module_name.split('.')[1]


    # Create log instance and get logger
    log = logging.getLogger()
    log.setLevel(logging.INFO)

    # Format logging output
    formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')


    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.setFormatter(formatter)

    # Set log directory in which all the logs will sit under.
    file_handler = logging.FileHandler(f"{log_dir}{module_name}_log.log")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    log.addHandler(file_handler)
    log.addHandler(stdout_handler)

