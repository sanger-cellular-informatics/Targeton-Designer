#!/usr/bin/env python3

import argparse
import sys

def get_args():
    parser = argparse.ArgumentParser(
        description = 'Targeton Designer CLI')

    parser.add_argument('command', type=str)

    return parser.parse_args()

def version_command():
    python_version = sys.version
    version = '0.0.1'

    print('Targeton Designer', 'version: ', version)
    print('Python version: ', python_version)

def resolve_command(command):
    if command == 'version':
        version_command()

def main():
    args = get_args()
    command = args.command

    resolve_command(command)


if __name__ == '__main__':
    main()