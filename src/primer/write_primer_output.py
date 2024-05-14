from collections import defaultdict
from typing import List

import pandas as pd
from os import path

from designer.output_data_classes import PrimerOutputData
from primer.designed_primer import DesignedPrimer
from primer.primer_pair import PrimerPair
from config.config import DesignerConfig
from utils.write_output_files import timestamped_dir, export_to_bed
from primer.filter.filter_response import PrimerPairDiscarded

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)

def write_primer_output(
    prefix='',
    primer_pairs=[],
    discarded_primer_pairs=[],
    existing_dir='',
    primer_type='LibAmp'
) -> PrimerOutputData:
    if existing_dir:
        export_dir = existing_dir
    else:
        export_dir = timestamped_dir(prefix)

    result = PrimerOutputData(export_dir)

    primer_rows = construct_primer_rows_bed_format(primer_pairs)
    result.bed = export_to_bed(primer_rows, export_dir)

    result.csv = export_primers_to_csv(primer_pairs, export_dir, primer_type)
    result.dir = export_dir

    logger.info(f"Primer files saved: {result.bed} {result.csv}")

    if discarded_primer_pairs:
        result.discarded_csv = export_discarded_primers_to_csv(
                                  discarded_primer_pairs,
                                  export_dir,
                                  primer_type)
        print(f"Discarded primer file saved: {result.discarded_csv}")
    else:
        print("No discarded primers")

    return result


def export_primers_to_csv(primer_pairs: List[PrimerPair], export_dir: str, primer_type: str) -> str:
    PRIMER3_OUTPUT_CSV = 'p3_output.csv'
    primers_csv_output_path = path.join(export_dir, PRIMER3_OUTPUT_CSV)

    primers_dataframe = _get_primers_dataframe(primer_pairs, primer_type)

    col_order = DesignerConfig().params['csv_column_order']
    primers_dataframe_ordered = _reorder_columns(col_order, primers_dataframe)
    primers_dataframe_ordered.to_csv(primers_csv_output_path, index=False)

    return primers_csv_output_path

def export_discarded_primers_to_csv(discarded_pairs: List[PrimerPairDiscarded],
                                    export_dir: str, primer_type: str) -> str:
    PRIMER3_DISCARDED_OUTPUT_CSV = 'discarded_pairs.csv'
    output_path = path.join(export_dir, PRIMER3_DISCARDED_OUTPUT_CSV)

    # create a data frame for output as csv
    discarded_df = _get_discarded_primer_dataframe(discarded_pairs, primer_type)
    col_order = DesignerConfig().params['csv_column_order']
    col_order.append('discard_reason')
    df_ordered = _reorder_columns(col_order, discarded_df)
    df_ordered.to_csv(output_path, index=False)
    return output_path

def _get_primers_dataframe(pairs: List[PrimerPair], primer_type: str) -> pd.DataFrame:
    primers_dict = defaultdict(list)

    for pair in pairs:

        for direction in ['forward', 'reverse']:
            primer = getattr(pair, direction)
            primers_dict['primer_type'].append(primer_type)
            primers_dict['primer'].append(primer.name)
            primers_dict['penalty'].append(primer.penalty)
            primers_dict['sequence'].append(primer.sequence)
            primers_dict['primer_start'].append(primer.primer_start)
            primers_dict['primer_end'].append(primer.primer_end)
            primers_dict['tm'].append(primer.tm)
            primers_dict['gc_percent'].append(primer.gc_percent)
            primers_dict['self_any_th'].append(primer.self_any_th)
            primers_dict['self_end_th'].append(primer.self_end_th)
            primers_dict['hairpin_th'].append(primer.hairpin_th)
            primers_dict['end_stability'].append(primer.end_stability)
        

        primers_dict['pair_uid'].extend([pair.uid] * 2)
        primers_dict['stringency'].extend([pair.stringency] * 2)
        primers_dict['chromosome'].extend([pair.chromosome] * 2)
        primers_dict['pre_targeton_start'].extend([pair.pre_targeton_start] * 2)
        primers_dict['pre_targeton_end'].extend([pair.pre_targeton_end] * 2)
        primers_dict['product_size'].extend([pair.product_size] * 2)
        primers_dict['targeton_id'].extend([pair.targeton_id] * 2)

    return pd.DataFrame(primers_dict).round(decimals=3)

def _get_discarded_primer_dataframe(discarded_pairs: List[PrimerPairDiscarded],
                                    primer_type: str) -> pd.DataFrame:
    discarded_primers_dict = defaultdict(list)
    for discarded_pair in discarded_pairs:
        primer_pair = discarded_pair.primer_pair
        _add_primer_pair(discarded_primers_dict, primer_pair, primer_type)
        discarded_primers_dict['discard_reason'].extend([discarded_pair.filter_applied]*2)

    return pd.DataFrame(discarded_primers_dict).round(decimals=3)

def _add_primer_pair(primers_dict: defaultdict(list),
                     pair: PrimerPair, primer_type: str) -> None:
    for direction in ['forward', 'reverse']:
        primer = getattr(pair, direction)
        primers_dict['primer_type'].append(primer_type)
        primers_dict['primer'].append(primer.name)
        primers_dict['penalty'].append(primer.penalty)
        primers_dict['sequence'].append(primer.sequence)
        primers_dict['primer_start'].append(primer.primer_start)
        primers_dict['primer_end'].append(primer.primer_end)
        primers_dict['tm'].append(primer.tm)
        primers_dict['gc_percent'].append(primer.gc_percent)
        primers_dict['self_any_th'].append(primer.self_any_th)
        primers_dict['self_end_th'].append(primer.self_end_th)
        primers_dict['hairpin_th'].append(primer.hairpin_th)
        primers_dict['end_stability'].append(primer.end_stability)

    primers_dict['pair_uid'].extend([pair.uid] * 2)
    primers_dict['stringency'].extend([pair.stringency] * 2)
    primers_dict['chromosome'].extend([pair.chromosome] * 2)
    primers_dict['pre_targeton_start'].extend([pair.pre_targeton_start] * 2)
    primers_dict['pre_targeton_end'].extend([pair.pre_targeton_end] * 2)
    primers_dict['product_size'].extend([pair.product_size] * 2)
    primers_dict['targeton_id'].extend([pair.targeton_id] * 2)

    return None

def _reorder_columns(csv_col_order: List[str],
                     dataframe: pd.DataFrame):
    
    col_order_unique = list(dict.fromkeys(csv_col_order))

    # raise Exception(col_order_unique)

    if not col_order_unique:
        logger.warning("Warning: empty csv_column_order list provided in config file, returning dataframe with default column order")
        return dataframe

    final_order = []
    for column in col_order_unique:
        if column not in dataframe.columns:
            logger.warning(f"Warning: '{column}' specified in config file not is not a column name")
        else:
            final_order.append(column)
    
    if not final_order:
        raise ValueError("All column names in config file are wrong")

    for column in dataframe.columns:
        if column not in final_order:
            logger.info(f"'{column}' column discarded as it is not in config file")

    return dataframe[final_order]


def construct_primer_rows_bed_format(pairs: List[PrimerPair]) -> list:
    primer_rows = []
    for pair in pairs:
        primer_rows.append(create_bed_row_for_primer(primer=pair.forward, chromosome=pair.chromosome))
        primer_rows.append(create_bed_row_for_primer(primer=pair.reverse, chromosome=pair.chromosome))

    return primer_rows


def create_bed_row_for_primer(primer: DesignedPrimer, chromosome: str) -> list:
    primer_row = [
        chromosome,
        primer.primer_start,
        primer.primer_end,
        primer.name,
        '0',
        primer.strand
    ]

    return primer_row
