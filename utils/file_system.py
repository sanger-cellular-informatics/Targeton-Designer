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


def write_to_text_file(dir_path, data, file_name):
    path = dir_path + '/' + file_name + '.txt'
   
    file_h = open(path, "wb")
    if isinstance(data, list):
        for row in data:
            file_h.write(row + "\n")
    else:
        file_h.write(data)
    file_h.close
    
    print('Wrote to file: ' + path)

    return path
