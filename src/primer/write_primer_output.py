from typing import List

import pandas as pd
from os import path

from designer.output_data_classes import PrimerOutputData
from primer.slice_data import SliceData
from primer.primer_pair import PrimerPair
from utils.write_output_files import timestamped_dir, export_to_bed


def write_primer_output(
    prefix='',
    primers=[],
    existing_dir='',
    primer_type='LibAmp'
) -> PrimerOutputData:

    if existing_dir:
        export_dir = existing_dir
    else:
        export_dir = timestamped_dir(prefix)

    result = PrimerOutputData(export_dir)

    bed_rows = construct_bed_format_from_pairs(primers)

    result.bed = export_to_bed(bed_rows, export_dir)
    result.csv = export_primers_to_csv(primers, export_dir, primer_type)
    result.dir = export_dir

    print('Primer files saved:', result.bed, result.csv)

    return result


def export_primers_to_csv(pairs: List[dict], export_dir: str, primer_type: str) -> str:
    PRIMER3_OUTPUT_CSV = 'p3_output.csv'

    #rows = construct_csv_format(slices)
    rows = construct_csv_format_from_pairs(pairs)
    if not rows.empty:
        rows.insert(0, 'primer_type', primer_type)

    full_path = path.join(export_dir, PRIMER3_OUTPUT_CSV)
    rows.to_csv(full_path, index=False)

    return full_path

def construct_csv_format_from_pairs(pairs: List[PrimerPair]) -> list:
    rows = pd.DataFrame()

    for pair in pairs:
        forward_df = transform_primer_to_df(pair.forward, pair.chromosome,
                                            pair.pre_targeton_start, pair.pre_targeton_end,
                                            pair.product_size, pair.targeton_id)
        reverse_df = transform_primer_to_df(pair.reverse, pair.chromosome,
                                            pair.pre_targeton_start, pair.pre_targeton_end,
                                            pair.product_size, pair.targeton_id)

        rows = pd.concat([rows, forward_df, reverse_df])

    return rows

def transform_primer_to_df(primer: dict, chromosome: str, pre_start: str, pre_end: str, prod_size: str,
                           targeton_id: str) -> pd.DataFrame:
    primer['chromosome'] = chromosome
    primer['pre_targeton_start'] = pre_start
    primer['pre_targeton_end'] = pre_end
    primer['product_size'] = prod_size
    primer['targeton_id'] = targeton_id

    primer.pop('coords', '')
    primer.pop('side', '')
    primer.pop('strand', '')
    primer.pop('pair_id', '')

    return pd.DataFrame([primer], [primer["primer"]])

def construct_bed_format_from_pairs(pairs: List[PrimerPair]) -> list:
    rows = []
    for pair in pairs:
        rows.append(create_bed_row_for_primer(pair.forward, pair.chromosome))
        rows.append(create_bed_row_for_primer(pair.reverse, pair.chromosome))
        
    return rows

def create_bed_row_for_primer(primer_data: dict, chromosome: str) -> dict:
    row = [
        chromosome,
        primer_data['primer_start'],
        primer_data['primer_end'],
        primer_data['primer'],
        '0',
        primer_data['strand']
    ]

    return row
