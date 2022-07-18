#!/usr/bin/env python3

import sys

from utils.arguments_parser import ParsedInputArguments
from utils.validate_files import validate_files
from utils.write_output_files import write_slicer_output, write_primer_output, write_ipcress_output
from slicer.slicer import Slicer
from primer.primer3 import Primer
from ipcress.ipcress import Ipcress


def version_command():
    python_version = sys.version
    version = '0.0.1'

    print('Targeton Designer version: ', version)
    print('Python version: ', python_version)


def slicer_command(args):
    validate_files(bed = args['bed'], fasta = args['fasta'])
    slicer = Slicer()
    slices = slicer.get_slices(args)

    return write_slicer_output(args['dir'], slices)


def primer_command(fasta, prefix = '', existing_dir = ''):
    validate_files(fasta = fasta)
    primer = Primer()
    primers = primer.get_primers(fasta)
    result = write_primer_output(
        primers = primers,
        prefix=prefix,
        existing_dir = existing_dir
    )

    return result


def ipcress_command(params, csv = '', existing_dir = ''):
    ipcress_input = params.copy()

    if csv:
        ipcress_input['p3_csv'] = csv
    if existing_dir:
        ipcress_input['dir'] = existing_dir

    validate_files(fasta = ipcress_input['fasta'],
                   txt = ipcress_input['primers'], p3_csv = ipcress_input['p3_csv'])

    ipcress = Ipcress()
    ipcress_result = ipcress.ipcress_runner(ipcress_input)

    write_ipcress_output(
        stnd = ipcress_result.stnd,
        err = ipcress_result.err,
        existing_dir = ipcress_input['dir']
    )


def resolve_command(args):
    command = args['command']

    if command == 'version':
        version_command()
    else:
        if command == 'slicer':
            slicer_command(args)

        if command == 'primer':
            primer_command(args['fasta'], prefix = args['dir'])

        if command == 'ipcress':
            ipcress_command(args)

        if command == 'design':
            slicer_result = slicer_command(args)
            primer_result = primer_command(slicer_result.fasta, existing_dir = slicer_result.dir)
            ipcress_command(args, csv = primer_result.csv,  existing_dir = slicer_result.dir)


def main():
    parsed_input = ParsedInputArguments()
    args = parsed_input.get_args()

    resolve_command(args)


if __name__ == '__main__':
    main()
