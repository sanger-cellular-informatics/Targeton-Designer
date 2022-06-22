import os

DEFAULT_FOLDER_NAME = 'output_folder'

class FolderCreator:
    dir = ''

    def set_dir(self, dir):
        self.dir = dir

    @classmethod
    def create(self):
        path = '/' + DEFAULT_FOLDER_NAME

        try:
            os.mkdir(path)
        except OSError as error:
            print(error)

        self.set_dir(self, path)