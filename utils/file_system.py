import os

DEFAULT_FOLDER_NAME = 'output_folder'

class FolderCreator:
    dir = ''

    def set_dir(self, dir):
        self.dir = dir

    @classmethod
    def get_dir(self):
        return(self.dir)

    @classmethod
    def create(self):
        path = './' + DEFAULT_FOLDER_NAME

        try:
            if not os.path.isdir(path):
                os.mkdir(path)
                print('Folder ' + path + ' is created')
            else:
                print('Folder ' + path + ' already exists')
        except OSError as error:
            print(error)

        self.set_dir(self, path)