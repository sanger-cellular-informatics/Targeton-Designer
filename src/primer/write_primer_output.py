from typing import List

import pandas as pd
from os import path

from designer.output_data_classes import PrimerOutputData
from primer.filter.designed_primer import DesignedPrimer
from primer.primer_pair import PrimerPair
from utils.write_output_files import timestamped_dir, export_to_bed


def write_primer_output(
    prefix='',
    primer_pairs=[],
    existing_dir='',
) -> PrimerOutputData:

    if existing_dir:
        export_dir = existing_dir
    else:
        export_dir = timestamped_dir(prefix)

    result = PrimerOutputData(export_dir)

    primer_rows = construct_primer_rows_bed_format(primer_pairs)
    result.bed = export_to_bed(primer_rows, export_dir)

    result.csv = export_primers_to_csv(primer_pairs, export_dir)
    result.dir = export_dir

    print('Primer files saved:', result.bed, result.csv)

    return result


def export_primers_to_csv(primer_pairs: List[PrimerPair], export_dir: str) -> str:
    PRIMER3_OUTPUT_CSV = 'p3_output.csv'
    primers_csv_output_path = path.join(export_dir, PRIMER3_OUTPUT_CSV)

    primers_dataframe = _get_primers_dataframe(primer_pairs)
    primers_dataframe.to_csv(primers_csv_output_path, index=False)

    return primers_csv_output_path

def _get_primers_dataframe(pairs: List[PrimerPair]) -> pd.DataFrame:
    primer_names = []
    primer_penalties = []
    primer_stringencies = []
    primer_sequences = []
    primer_starts = []
    primer_ends = []
    primer_temperatures_melting = []
    primer_gc_percentages = []
    primer_self_any_ths = []
    primer_self_end_ths = []
    primer_hairpin_ths = []
    primer_end_stabilities = []
    pair_chromosomes = []
    pair_pre_targeton_starts = []
    pair_pre_targeton_ends = []

    for pair in pairs:
        forward_primer = pair.forward
        reserve_primer = pair.reverse

        primer_names.extend([forward_primer.name, reserve_primer.name])
        primer_penalties.extend([forward_primer.penalty, reserve_primer.penalty])
        primer_stringencies.extend([forward_primer.stringency, reserve_primer.stringency])
        primer_sequences.extend([forward_primer.sequence, reserve_primer.sequence])
        primer_starts.extend([forward_primer.primer_start, reserve_primer.primer_start])
        primer_ends.extend([forward_primer.primer_end, reserve_primer.primer_end])
        primer_temperatures_melting.extend([forward_primer.tm, reserve_primer.tm])
        primer_gc_percentages.extend([forward_primer.gc_percent, reserve_primer.gc_percent])
        primer_self_any_ths.extend([forward_primer.self_any_th, reserve_primer.self_any_th])
        primer_self_end_ths.extend([forward_primer.self_end_th, reserve_primer.self_end_th])
        primer_hairpin_ths.extend([forward_primer.hairpin_th, reserve_primer.hairpin_th])
        primer_end_stabilities.extend([forward_primer.end_stability, reserve_primer.end_stability])

        pair_chromosomes.extend([pair.chromosome] * 2)
        pair_pre_targeton_starts.extend([pair.pre_targeton_start] * 2)
        pair_pre_targeton_ends.extend([pair.pre_targeton_end] * 2)

    all_primers_data = {
        'primer': primer_names,
        'penalty': primer_penalties,
        'stringency': primer_stringencies,
        'sequence': primer_sequences,
        'primer_start': primer_starts,
        'primer_end': primer_ends,
        'tm': primer_temperatures_melting,
        'gc_percent': primer_gc_percentages,
        'self_any_th': primer_self_any_ths,
        'self_end_th': primer_self_end_ths,
        'hairpin_th': primer_hairpin_ths,
        'end_stability': primer_end_stabilities,
        'chromosome': pair_chromosomes,
        'pre_targeton_start': pair_pre_targeton_starts,
        'pre_targeton_end': pair_pre_targeton_ends
    }

    return pd.DataFrame(all_primers_data)


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
