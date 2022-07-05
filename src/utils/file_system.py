import os
from datetime import datetime

class FolderCreator:
    dir = ''

    def set_dir(self, dir):
        self.dir = dir

    @classmethod
    def get_dir(self):
        return(self.dir)

    @classmethod
    def create(self, path):
        try:
            if not os.path.isdir(path):
                os.mkdir(path)
                print(f'Folder {path} is created')
            else:
                print(f'Warning: {path} already exists and files may be overwritten')
        except OSError as error:
            raise FolderCreatorError(f'Unexpected OSError: {error}')

        self.set_dir(self, path)

    @classmethod
    def create_timestamped(self, prefix):
        self.create(f'{prefix}_' + datetime.now().strftime('%Y%m%d%H%M%S%f'))


def check_file_exists(file):
    if not os.path.exists(file):
        raise FileNotFoundError(f'Unable to find file: {file}')
