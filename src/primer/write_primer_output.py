from collections import defaultdict
from dataclasses import dataclass
from typing import List, Optional

import pandas as pd
from os import path

from config.ipcress_params import IpcressParameters
from designer.output_data_classes import PrimerOutputData
from primer.designed_primer import DesignedPrimer
from primer.primer_pair import PrimerPair
from utils.write_output_files import timestamped_dir, export_to_bed
from primer.filter.filter_response import PrimerPairDiscarded

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)

# Define constants
IPCRESS_COLUMN_ORDER = ["id", "forward_sequence", "reverse_sequence", "min_size", "max_size"]
PRIMER3_OUTPUT_CSV = 'p3_output.csv'
PRIMER3_DISCARDED_OUTPUT_CSV = 'discarded_pairs.csv'
OPTIMAL_PRIMERS_CSV = 'optimal_primer_pairs.csv'
PRIMER_PAIRS_FOR_IPCRESS_CSV = 'primer_pairs_for_ipcress.csv'


def write_primer_output(
    sorted_primer_pairs: List[PrimerPair],
    column_order: List[str],
    prefix='',
    primer_pairs=[],
    discarded_primer_pairs=[],
    existing_dir='',
    primer_type='LibAmp',
    ipcress_params: Optional[IpcressParameters] = None,
) -> PrimerOutputData:
    export_dir = existing_dir or timestamped_dir(prefix)
    result = PrimerOutputData(export_dir)

    # DataFrame of primers, each primer on its own row
    primers_df = _get_primers_dataframe(sorted_primer_pairs, primer_type)

    if primer_pairs:
        primer_rows = construct_primer_rows_bed_format(primer_pairs)
        result.bed = export_to_bed(primer_rows, export_dir)

        result.csv = export_primers_to_csv(primers_df, export_dir, column_order)
        result.optimal_primer_pairs_csv = export_three_optimal_primer_pairs_to_csv(
            primers_df,
            export_dir,
            column_order
        )

        logger.info(f"Primer files saved: {result.bed}, {result.csv}, {result.optimal_primer_pairs_csv}")

        if ipcress_params and ipcress_params.write_ipcress_file:

            # DataFrame of primer pairs, forward and reverse primer on the same row
            primer_pairs_df = _get_primer_pairs_dataframe_for_ipcress(
                sorted_primer_pairs,  
                ipcress_params.min_size,
                ipcress_params.max_size,
            )
            
            result.primer_pairs_for_ipcress = export_pairs_for_ipcress_to_csv(
                primer_pairs_df, 
                export_dir, 
                IPCRESS_COLUMN_ORDER
            )
        
            logger.info(f"Primer pairs for ipcress saved: {result.primer_pairs_for_ipcress}")


    if discarded_primer_pairs:
        discarded_primers_df = _get_discarded_primer_dataframe(discarded_primer_pairs, primer_type)
        column_order.append('discard_reason')

        result.discarded_csv = export_discarded_primers_to_csv(
                                  discarded_primers_df,
                                  export_dir,
                                  column_order
                                  )
        logger.info(f"Discarded primer file saved: {result.discarded_csv}")
    else:
        logger.info("No discarded primers")

    return result
    


def export_primers_to_csv(primers_dataframe: pd.DataFrame, export_dir: str, column_order: List[str]) -> str:
    primers_csv_output_path = path.join(export_dir, PRIMER3_OUTPUT_CSV)

    write_dataframe_to_csv(primers_dataframe, column_order, primers_csv_output_path)

    return primers_csv_output_path

def export_discarded_primers_to_csv(discarded_df: pd.DataFrame, export_dir: str, column_order: List[str]) -> str:
    output_path = path.join(export_dir, PRIMER3_DISCARDED_OUTPUT_CSV)

    write_dataframe_to_csv(discarded_df, column_order, output_path)

    return output_path

def export_three_optimal_primer_pairs_to_csv(df: pd.DataFrame, export_dir: str, column_order: List[str]) -> str:
    primers_csv_output_path = path.join(export_dir, OPTIMAL_PRIMERS_CSV)

    optimal_primers_df = df.head(6)
    
    if len(optimal_primers_df.index) < 6:
        logger.warning("Less than 3 primer pairs returned by Primer3")

    write_dataframe_to_csv(optimal_primers_df, column_order, primers_csv_output_path)

    return primers_csv_output_path

def export_pairs_for_ipcress_to_csv(df: pd.DataFrame, export_dir: str, column_order: List[str]) -> str:
    primers_csv_output_path = path.join(export_dir, PRIMER_PAIRS_FOR_IPCRESS_CSV)

    write_dataframe_to_csv(
        df,
        column_order,
        primers_csv_output_path,
        include_header=False,
        sep="\t",
    )

    return primers_csv_output_path

def write_dataframe_to_csv(
    df: pd.DataFrame,
    cols: List[str],
    output_path: str,
    include_header: bool = True,
    sep: str = ",",
) -> None:
    df_ordered = _reorder_columns(cols, df)
    df_ordered.to_csv(output_path, index=False, header=include_header, sep=sep)
    return None

def _get_primers_dataframe(pairs: List[PrimerPair], primer_type: str) -> pd.DataFrame:
    primers_dict = defaultdict(list)

    for pair in pairs:
        _add_primer_pair(primers_dict, pair, primer_type)

    return pd.DataFrame(primers_dict).round(decimals=3)

def _get_discarded_primer_dataframe(discarded_pairs: List[PrimerPairDiscarded],
                                    primer_type: str) -> pd.DataFrame:
    discarded_primers_dict = defaultdict(list)
    for discarded_pair in discarded_pairs:
        _add_primer_pair(discarded_primers_dict, discarded_pair, primer_type)
        discarded_primers_dict['discard_reason'].extend([discarded_pair.reason_discarded]*2)

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

def _get_primer_pairs_dataframe_for_ipcress(
        pairs: List[PrimerPair], 
        min_size: str, 
        max_size: str,
    ) -> pd.DataFrame:

    primer_pairs_dict = defaultdict(list)

    for pair in pairs:
        primer_pairs_dict['id'].append(_normalize_primer_name_to_pair_id(pair.forward.name))
        primer_pairs_dict['forward_sequence'].append(pair.forward.sequence)
        primer_pairs_dict['reverse_sequence'].append(pair.reverse.sequence)
        primer_pairs_dict['min_size'].append(min_size)
        primer_pairs_dict['max_size'].append(max_size)

    return pd.DataFrame(primer_pairs_dict).round(decimals=3)

def _normalize_primer_name_to_pair_id(primer_name: str) -> str:
    match = primer_name.rsplit("_", 1)
    if len(match) != 2:
        return primer_name

    head, pair_number = match
    if head.endswith("F") or head.endswith("R"):
        head = head[:-1]
    return f"{head}_{pair_number}"

def _reorder_columns(csv_col_order: List[str],
                     dataframe: pd.DataFrame):

    col_order_unique = list(dict.fromkeys(csv_col_order))
    _check_unique_columns(col_order_unique, csv_col_order)

    if not col_order_unique:
        logger.warning("Empty csv_column_order list provided in config file, returning dataframe with default column order")
        return dataframe

    final_order = []
    for column in col_order_unique:
        if column not in dataframe.columns:
            logger.warning(f"'{column}' specified in config file not is not a column name")
        else:
            final_order.append(column)

    if not final_order:
        raise ValueError("All column names in config file are wrong")

    for column in dataframe.columns:
        if column not in final_order:
            logger.warning(f"'{column}' column discarded as it is not in config file")

    return dataframe[final_order]


def _check_unique_columns(col_order_unique, csv_col_order):
    if len(col_order_unique) != len(csv_col_order):
        discarded = []
    for column in csv_col_order:
        if csv_col_order.count(column) > 1 and column not in discarded:
            discarded.append(column)
            logger.warning(f"'{column}' duplicated in config file, only first instance retained")


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
