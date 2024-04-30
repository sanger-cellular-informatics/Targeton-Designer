import logging
from logging import Logger
import sys, os


class CustomLogger(Logger):

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.name = name
    


    def logs_to_log(self) -> None:

        """
        Function logs_to_log() will log out all of current logging output to a log file,
        which sits under "/src/logs" directory. It will grab the logs outputted from a file
        in which it's instantiated is added.

        This function take one argument module_name and returns None.  
        For example: "ParentModule.Module". 

        Here, default logger level is DEBUG. 
        This means you can log INFO, ERROR, WARNING and EXCEPTION.
        """


        module_name = self.name.split('.')[1]
        
        # Get directory for logs
        log_directory = os.path.join(os.getcwd(),"src/logs/")

        if not os.path.exists(log_directory):
            os.makedirs(log_directory)
        log_file_path = os.path.join(log_directory, f"{module_name}_log.log")
    
        # Create log instance and get logger
        log = logging.getLogger()
        log.setLevel(logging.INFO)

        # Format logging output
        logging.basicConfig(
                            filename=log_file_path, 
                            filemode='a', 
                            level=logging.DEBUG,
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                            force=True
        )

        return '%(asctime)s | %(name)s | %(levelname)s | %(message)s'

