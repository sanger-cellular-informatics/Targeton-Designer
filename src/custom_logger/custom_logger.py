import coloredlogs
import logging
import os


def setup_logger(module_name: str) -> logging.Logger:
    """
    
        Function setup_logger() will log out all of current logging output to a log file,
        which sits under "/src/logs/logs.log" directory. It will grab the logs outputted from a 
        file in which it's instantiated is added.

        This function take one argument module_name and returns Logger instance.  
        For example: setup_logger(ParentModule.Module: str) -> logging.Logger

        Here, default logger level is DEBUG. 
        This means you can log INFO, WARNING, ERROR, and EXCEPTION.

        Example of use cases:

        import setup_logger

        logger = setup_logger(__name__)
        
        logger.info("Info")
        logger.warn("Warning!")
        logger.error("Error")
        logger.exception("Exception")

    """

    # Create a logger object.
    logger = logging.getLogger(module_name)
        
    # Get directory for logs
    log_directory = os.path.join(os.getcwd(), "src/logs/")

    # Ensure the log directory exists
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)
    log_file_path = os.path.join(log_directory, f"logs.log")

    # Set up logging to file and format
    logging.basicConfig(
        filename=log_file_path, 
        filemode='a', 
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        force=True,
    )

    # Install colored logs for console output
    coloredlogs.install(level='DEBUG', logger=logger)

    # Return configured logging instance
    return logger

