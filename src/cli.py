#!/usr/bin/env python3

import sys

from utils.arguments_parser import ParsedInputArguments
from utils.validate_files import validate_files
from utils.write_output_files import write_slicer_output, write_primer_output
from slicer.slicer import main as slicer
from primer.primer3 import main as primer

def version_command():
    python_version = sys.version
    version = '0.0.1'

    print('Targeton Designer', 'version: ', version)
    print('Python version: ', python_version)

def slicer_command(args):
    validate_files(bed = args['bed'], fasta = args['fasta'])
    slices = slicer(args)
    created_output_dir = write_slicer_output(args['dir'], slices)

    return created_output_dir

def primer_command(args, output_dir = ''):
    validate_files(fasta = args['fasta'])
    primers = primer(args['fasta'])
    write_primer_output(prefix = args['dir'], primers = primers, existing_dir = output_dir)

def resolve_command(args):
    command = args['command']

    if command == 'version':
        version_command()
    else:
        print('args: ', args)

        if command == 'slicer':
            slicer_command(args)

        if command == 'primer':
            primer_command(args)

        if command == 'design':
            output_dir = slicer_command(args)
            primer_command(args, output_dir)

def main():
    parsed_input = ParsedInputArguments()
    args = parsed_input.get_args()

    resolve_command(args)

if __name__ == '__main__':
    main()