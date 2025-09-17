#!/usr/bin/env python3
import sys
from os import path

from config.config import DesignerConfig
from primer.filter.filter_manager import FilterManager
from primer.slice_data import SliceData
from utils.arguments_parser import ParsedInputArguments
from utils.validate_files import validate_files, validate_fasta_format
from utils.write_output_files import (
    write_slicer_output,
    write_targeton_csv,
    write_scoring_output,
    write_primer_design_output,
)

from designer.output_data_classes import (
    SlicerOutputData,
    PrimerOutputData,
    ScoringOutputData,
    PrimerDesignerOutputData,
    DesignOutputData,
)
from primer.write_primer_output import write_primer_output
from slicer.slicer import Slicer
from primer.primer3 import Primer3
from primer_designer import PrimerDesigner
from post_primer_pairs import post_primer_pairs
from primer.ranker.ranker import Ranker

sys.path.append(path.abspath(path.join(path.dirname(__file__), '../sge-primer-scoring/src')))
from scoring import Scoring


def version_command():
    python_version = sys.version
    version = '0.0.1'

    print('Primer Designer version: ', version)
    print('Python version: ', python_version)


def slicer_command(args) -> SlicerOutputData:
    validate_files(bed=args['bed'], fasta=args['fasta'])
    slicer = Slicer()
    slices = slicer.get_slices(args)

    return write_slicer_output(args['dir'], slices)


def primer_command(
        args: dict
) -> PrimerOutputData:
    PRIMER_TYPE = 'LibAmp'
    config = DesignerConfig(args)

    validate_fasta_format(config.fasta)

    slice_data = SliceData.get_first_slice_data(config.fasta, config.region_padding, config.region_avoid)

    primers = Primer3(config.stringency_vector, config.primer3_params).get_primers(slice_data)

    filters_response = FilterManager(config.filters).apply_filters(primers)

    ranked_primer_pairs_df = (Ranker(config.ranking)
                              .rank(primer_type=PRIMER_TYPE, primer_pairs=filters_response.primer_pairs_to_keep))

    primer_result = write_primer_output(
        primer_pairs_df=ranked_primer_pairs_df,
        primer_pairs=filters_response.primer_pairs_to_keep,
        discarded_primer_pairs=filters_response.primer_pairs_to_discard,
        prefix=config.prefix_output_dir,
        primer_type=PRIMER_TYPE,
        column_order=config.csv_column_order
    )

    return primer_result


def collate_primer_designer_data_command(
    design_output_data : DesignOutputData,
    primer_designer=PrimerDesigner(),
    prefix='',
    existing_dir=''
) -> PrimerDesignerOutputData:
    validate_files(p3_csv=design_output_data.p3_csv, score_tsv=design_output_data.scoring_tsv)

    primer_designer.from_design_output(design_output_data)
    primer_designer_result = write_primer_design_output(
        primer_designer,
        prefix=prefix,
        existing_dir=existing_dir,
    )
    return primer_designer_result


def scoring_command(ipcress_output, mismatch, output_tsv, targeton_csv=None) -> ScoringOutputData:
    scoring = Scoring(ipcress_output, mismatch, targeton_csv)
    scoring.add_scores_to_df()

    result = write_scoring_output(scoring, output_tsv)

    return result


def design_command(args) -> DesignOutputData:
    primer_result = primer_command(args=args)

    output_dir = primer_result.dir

    design_result = DesignOutputData(output_dir)
    # Primer
    design_result.p3_bed = primer_result.bed
    design_result.p3_csv = primer_result.csv

    return design_result


def post_primers(primer_json) -> None:
    validate_files(primer_json=primer_json)
    post_primer_pairs(primer_json)


def resolve_command(args):
    command = args['command']

    if command == 'version':
        version_command()
    else:
        if command == 'slicer':
            slicer_command(args)

        if command == 'primer':
            primer_command(args=args)

        if command == 'collate_primer_data':
            design_output_data = DesignOutputData()
            design_output_data.p3_csv = args['p3_csv']
            design_output_data.scoring_tsv = args['score_tsv']
            collate_primer_designer_data_command(design_output_data, prefix=args['dir'])

        if command == 'generate_targeton_csv':
            write_targeton_csv(args['primers'], args['bed'], args['dir'])

        if command == 'scoring':
            scoring_command(
                args['ipcress_file'],
                args['scoring_mismatch'],
                args['output_tsv'],
                args['targeton_csv'],
            )

        if command == 'design':
            design_command(args)

        if command == 'post_primers':
            post_primers(args['primer_json'])


def main():
    parsed_input = ParsedInputArguments()
    args = parsed_input.get_args()

    resolve_command(args)


if __name__ == '__main__':
    main()
