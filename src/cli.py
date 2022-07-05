#!/usr/bin/env python3

import sys

from utils.arguments_parser import ParsedInputArguments
from utils.validate_files import validate_files
from utils.write_output_files import write_slicer_output
from slicer.slicer import main as slicer

def version_command():
    python_version = sys.version
    version = '0.0.1'

    print('Targeton Designer', 'version: ', version)
    print('Python version: ', python_version)

def slicer_command(args):
    validate_files(args['bed'], args['fasta'])
    slices = slicer(args)
    write_slicer_output(args['dir'], slices)

def resolve_command(args):
    command = args['command']

    if command == 'version':
        version_command()
    else:
        print('args: ', args)

        if command == 'slicer':
            slicer_command(args)

def main():
    parsed_input = ParsedInputArguments()
    args = parsed_input.get_args()

    resolve_command(args)

if __name__ == '__main__':
    main()