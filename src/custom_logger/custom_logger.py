# import logging
# from logging import Logger
# import os
import coloredlogs
import logging
import os


def setup_logger(module_name: str) -> logging.Logger:
    """
        Function setup_logger() will log out all of current logging output to a log file, \n
        which sits under "/src/logs" directory. 
        It will grab the logs outputted from a file in which it's instantiated is added.

        This function take one argument module_name and returns Logger instance.  
        For example: (ParentModule.Module: str) -> logging.Logger

        Here, default logger level is DEBUG. 
        This means you can log INFO, ERROR, WARNING and EXCEPTION.
    """
    
    # Create a logger object.
    logger = logging.getLogger(module_name)

    # Determine the base module name if used as a package, otherwise just use the name
    module_name = module_name.split('.')[1] if '.' in module_name else module_name
        
    # Get directory for logs
    log_directory = os.path.join(os.getcwd(), "src/logs/")

    # Ensure the log directory exists
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)
    log_file_path = os.path.join(log_directory, f"{module_name}_log.log")

    # Set up logging to file and format
    logging.basicConfig(
        filename=log_file_path, 
        filemode='a', 
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        force=True
    )

    # Install colored logs for console output
    coloredlogs.install(level='DEBUG', logger=logger)

    return logger

