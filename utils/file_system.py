import csv

from os import path, mkdir
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
            if not path.isdir(path):
                mkdir(path)
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
    if not path.exists(file):
        raise FileNotFoundError(f'Unable to find file: {file}')


def write_to_text_file(dir_path, data, file_name):
    path = dir_path + '/' + file_name + '.txt'
   
    file_h = {}
    if isinstance(data, list):
        file_h = open(path, "w")
        for row in data:
            file_h.write(row + "\n")
    else:
        file_h = open(path, "wb")
        file_h.write(data)
    file_h.close
    
    print('Wrote to file: ' + path)

    return path


def read_csv_to_dict(csv_path):
    check_file_exists(csv_path)
    
    data = []
    with open(csv_path) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')
        for row in reader:
            data.append(row)
     
    return data
