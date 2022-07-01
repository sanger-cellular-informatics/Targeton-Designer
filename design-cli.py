#!/usr/bin/env python3

import argparse
import sys

def get_args():
    parser = argparse.ArgumentParser(
        description = 'Targeton Designer CLI')

    parser.add_argument('command', type=str)

    return parser.parse_args()

def version_command():
    td_version = '0.0.0'
    python_version = sys.version

    print('Targeton Designer', 'version: ', td_version)
    print('Python version: ', python_version)

def main():
    args = get_args()
    command = args.command

    if command == 'version':
        version_command()

if __name__ == '__main__':
    main()