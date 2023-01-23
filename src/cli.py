#!/usr/bin/env python3

import sys
import os

from utils.arguments_parser import ParsedInputArguments
from utils.validate_files import validate_files
from utils.write_output_files import write_slicer_output, write_primer_output, write_ipcress_input, write_ipcress_output
from slicer.slicer import Slicer
from primer.primer3 import Primer
from ipcress.ipcress import Ipcress
from adapters.primer3_to_ipcress import Primer3ToIpcressAdapter

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../sge-primer-scoring/src'))
)
from scoring import Scoring


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


def primer_command(fasta = '', prefix = '', existing_dir = ''):
    validate_files(fasta = fasta)
    primer = Primer()
    primers = primer.get_primers(fasta)

    result = write_primer_output(
        primers = primers,
        prefix = prefix,
        existing_dir = existing_dir
    )

    return result


def primer_for_ipcress(fasta = '', prefix = '', min = 0, max = 0):
    primer_result = primer_command(fasta = fasta, prefix = prefix)

    adapter = Primer3ToIpcressAdapter()
    adapter.prepare_input(primer_result.csv, min, max, primer_result.dir)

    input_path = write_ipcress_input(primer_result.dir, adapter.formatted_primers)

    result = adapter.formatted_primers

    return result


def ipcress_command(params, csv = '', existing_dir = ''):
    ipcress_params = params.copy()

    if csv:
        ipcress_params['p3_csv'] = csv
    if existing_dir:
        ipcress_params['dir'] = existing_dir

    validate_files(
        fasta = ipcress_params['fasta'],
        txt = ipcress_params['primers'],
        p3_csv = ipcress_params['p3_csv']
    )

    ipcress = Ipcress(ipcress_params)
    ipcress_result = ipcress.run()

    write_ipcress_output(
        stnd = ipcress_result.stnd,
        err = ipcress_result.err,
        existing_dir = ipcress_params['dir']
    )


def scoring_command(ipcress_output, mismatch, output_tsv, targeton_csv):
    scoring = Scoring(ipcress_output, mismatch, targeton_csv)
    scoring.add_scores_to_df()
    scoring.save_mismatches(output_tsv)


def resolve_command(args):
    command = args['command']

    if command == 'version':
        version_command()
    else:
        if command == 'slicer':
            slicer_command(args)

        if command == 'primer':
            primer_command(fasta = args['fasta'], prefix = args['dir'])

        if command == 'primer_for_ipcress':
            primer_for_ipcress(fasta = args['fasta'], prefix = args['dir'], min = args['min'], max = args['max'])

        if command == 'ipcress':
            ipcress_command(args)

        if command == 'scoring':
            scoring_command(
                args['ipcress_file'],
                args['scoring_mismatch'],
                args['output_tsv'],
                args['targeton_csv'],
            )

        if command == 'design':
            slicer_result = slicer_command(args)
            primer_result = primer_command(fasta = slicer_result.fasta, existing_dir = slicer_result.dir)
            ipcress_command(args, csv = primer_result.csv,  existing_dir = slicer_result.dir)


def main():
    parsed_input = ParsedInputArguments()
    args = parsed_input.get_args()

    resolve_command(args)


if __name__ == '__main__':
    main()
