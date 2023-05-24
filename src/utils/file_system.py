import csv
import json
from os import path, makedirs
from datetime import datetime
from typing import List, Any
import re

from utils.exceptions import FolderCreatorError, FileFormatError


class FolderCreator:
    dir = ''

    def set_dir(self, dir: str) -> None:
        self.dir = dir

    @classmethod
    def get_dir(self) -> str:
        return (self.dir)

    @classmethod
    def create(self, dir: str) -> None:
        try:
            if not path.isdir(dir):
                makedirs(dir)
                print(f'Folder {dir} is created')
            else:
                print(f'Warning: {dir} already exists and files may be overwritten')
        except OSError as error:
            raise FolderCreatorError(f'Unexpected OSError: {error}')

        self.set_dir(self, dir)

    @classmethod
    def create_timestamped(self, parent: str) -> None:
        PREFIX = 'td'

        dir_name = f'{PREFIX}_' + datetime.now().strftime('%Y%m%d%H%M%S%f')
        dir_full_path = path.join(parent, dir_name)

        self.create(dir_full_path)


def check_file_exists(file):
    if not path.exists(file):
        raise FileNotFoundError(f'Unable to find file: {file}')


def write_to_text_file(dir_path, data, file_name):
    file_path = dir_path + '/' + file_name + '.txt'

    file_h = {}
    if isinstance(data, list):
        file_h = open(file_path, "w")
        for row in data:
            file_h.write(row + "\n")
    else:
        file_h = open(file_path, "wb")
        file_h.write(data)
    file_h.close()

    print('Wrote to file: ' + file_path)

    return file_path


def read_csv_to_list_dict(csv_path, delimiter=',') -> List[dict]:
    check_file_exists(csv_path)

    data = []
    with open(csv_path, newline='', mode='r') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=delimiter)
        for row in reader:
            for key, value in row.items():
                row[key] = convert_string_to_python(value)
            data.append(row)

    return data


def read_csv(csv_path, delimiter=',') -> Any:
    check_file_exists(csv_path)
    data = []
    with open(csv_path, mode='r', newline='') as csv_file:
        reader = csv.reader(csv_file, delimiter=delimiter)
        for row in reader:
            new_row = []
            for element in row:
                new_row.append(convert_string_to_python(element))
            data.append(new_row)
    return data


def convert_string_to_python(string: str) -> Any:
    # group 0: all, 1: first element if iterable, 2: inner element no delimiter
    patterns = {
        list: r"^\[.+]$",  # list start/end with []
        tuple: r"^\(.+\)$",  # tuple start/end with ()
        dict: r"^\{.+\}$",  # dict start/end with {} with : somewhere between
        int: r"\d+",  # Integer
        float: r"[0-9.]+",  # Numeric only
        bool: r"True|False"  # Bool
    }
    for dtype, pattern in patterns.items():
        match = re.search(pattern, string)
        if match:
            if dtype in (list, tuple):
                split_string = string.split(',')
                data = []
                for element in split_string:
                    data.append(convert_string_to_python(element))
            elif dtype == dict:
                split_string = string.split(',')
                data = dict()
                for element in split_string:
                    key, value = element.split(':')
                    data.update({key: convert_string_to_python(value)})
            else:
                data = match.group()
            return dtype(data)
    return string  # default string


def regex_findall(pattern: str, string: str) -> str:
    matches = re.findall(pattern, string)
    matches = [match for match in matches]
    return matches


def parse_json(file_path: str) -> dict:
    with open(file_path, "r") as file:
        try:
            result = json.load(file)
        except Exception as err:
            raise FileFormatError

    return result
