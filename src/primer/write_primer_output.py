from typing import List

import pandas as pd
from os import path

from designer.output_data_classes import PrimerOutputData
from primer.slice_data import SliceData
from primer.primer_pair import PrimerPair
from config.config import DesignerConfig
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
    rows_unordered = construct_csv_format_from_pairs(pairs)
    if not rows_unordered.empty:
        rows_unordered.insert(0, 'primer_type', primer_type)

        col_order = DesignerConfig().params['csv_column_order']

        # col_order_default will have to be changed to a default config file stored somewhere
        col_order_default = ['primer_type',
                            'primer',
                            'penalty', 
                            'stringency', 
                            'sequence',
                            'primer_start',
                            'primer_end',
                            'tm',
                            'gc_percent',
                            'self_any_th',
                            'self_end_th',
                            'hairpin_th',
                            'end_stability',
                            'chromosome',
                            'pre_targeton_start',
                            'pre_targeton_end',
                            'product_size']    
        col_order = check_column_order_config(col_order, col_order_default)

        rows = rows_unordered[col_order]

    full_path = path.join(export_dir, PRIMER3_OUTPUT_CSV)
    rows.to_csv(full_path, index=False)

    return full_path

def construct_csv_format_from_pairs(pairs: List[PrimerPair]) -> list:
    rows = pd.DataFrame()

    for pair in pairs:
        forward_df = transform_primer_to_df(pair.forward, pair.chromosome,
                                            pair.pre_targeton_start, pair.pre_targeton_end,
                                            pair.product_size)
        reverse_df = transform_primer_to_df(pair.reverse, pair.chromosome,
                                            pair.pre_targeton_start, pair.pre_targeton_end,
                                            pair.product_size)

        rows = pd.concat([rows, forward_df, reverse_df])

    return rows

def transform_primer_to_df(primer: dict, chromosome: str, pre_start: str, pre_end: str, prod_size: str) -> pd.DataFrame:
    primer['chromosome'] = chromosome
    primer['pre_targeton_start'] = pre_start
    primer['pre_targeton_end'] = pre_end
    primer['product_size'] = prod_size

    primer.pop('coords', '')
    primer.pop('side', '')
    primer.pop('strand', '')
    primer.pop('pair_id', '')

    return pd.DataFrame([primer], [primer["primer"]])

def check_column_order_config(
        col_order: list,
        col_order_default: list
) -> list:
    # Empty
    if col_order == []:
        print("Error: csv_column_order list empty in config file") # Will have to change to logging

    # All incorrect columns
    elif (not any(x in col_order for x in col_order_default)):
        print("Error: no csv_column_order elements in config file match csv column names")

    else:
        # Some incorrect columns
        if len(set(col_order) - set(col_order_default)) > 0:
            print("Warning: ignored csv_column_order elements in config file that don't match csv column names")
            col_order = [x for x in col_order if x in col_order_default]

        # Columns are missing
        if not set(col_order_default).issubset(col_order):
            print("Warning: missing columns in config file were not used in csv output")
        
        # Duplicates
        if len(col_order) != len(set(col_order)):
            print("Warning: only first intance of duplicate csv_column_order elements from config file retained")
            col_order = list(dict.fromkeys(col_order))
    
    return col_order

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
