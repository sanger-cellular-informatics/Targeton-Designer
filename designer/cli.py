#!/usr/bin/env python3

import sys
from arguments_parser import ParsedInputArguments
from slicer import main as run_slicer

def version_command():
    python_version = sys.version
    version = '0.0.1'

    print('Targeton Designer', 'version: ', version)
    print('Python version: ', python_version)

def resolve_command(args):
    command = args['command']

    if command == 'version':
        version_command()
    else:
        print('args: ', args)

        if command == 'slicer':
            run_slicer(args)

def main():
    parsed_input = ParsedInputArguments()
    args = parsed_input.get_args()

    resolve_command(args)

if __name__ == '__main__':
    main()