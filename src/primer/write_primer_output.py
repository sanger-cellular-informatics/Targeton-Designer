from typing import List

import pandas as pd
from os import path

from designer.output_data_classes import PrimerOutputData
from primer.slice_data import SliceData
from utils.write_output_files import timestamped_dir, export_to_bed


def write_primer_output(
    prefix='',
    primers=[],
    existing_dir='',
) -> PrimerOutputData:

    if existing_dir:
        export_dir = existing_dir
    else:
        export_dir = timestamped_dir(prefix)

    result = PrimerOutputData(export_dir)

    bed_rows = construct_bed_format(primers)

    result.bed = export_to_bed(bed_rows, export_dir)
    result.csv = export_primers_to_csv(primers, export_dir)
    result.dir = export_dir

    print('Primer files saved:', result.bed, result.csv)

    return result


def export_primers_to_csv(slices: List[dict], export_dir: str) -> str:
    PRIMER3_OUTPUT_CSV = 'p3_output.csv'

    rows = construct_csv_format(slices)

    full_path = path.join(export_dir, PRIMER3_OUTPUT_CSV)
    rows.to_csv(full_path, index=False)

    return full_path


def construct_csv_format(slices: List[SliceData]) -> list:
    rows = pd.DataFrame()

    for slice_data in slices:
        primers = slice_data.primers
        for primer in primers:
            primers[primer]['chr'] = slice_data.chrom

            primers[primer].pop('coords', '')
            primers[primer].pop('side', '')
            primers[primer].pop('strand', '')
            primers[primer].pop('pair_id', '')

            primer_df = pd.DataFrame(primers[primer], [primer])
            rows = pd.concat([rows, primer_df])

    return rows.drop_duplicates()


def construct_bed_format(slices: List[SliceData]) -> list:
    rows = []
    for slice_data in slices:
        primers = slice_data.primers
        for primer in primers:
            primer_data = primers[primer]
            # chr,chrStart,chrEnd,name,score,strand
            # Score unknown until iPCRess
            row = [slice_data.chrom, primer_data['primer_start'],
                primer_data['primer_end'], primer, '0', primer_data['strand']]
            rows.append(row)
    return rows