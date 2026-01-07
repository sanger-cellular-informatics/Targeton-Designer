from __future__ import annotations

import csv
import re

from typing import TYPE_CHECKING, List, Union
from os import path
import os
from pathlib import Path

from pybedtools import BedTool
from utils.file_system import FolderCreator
from utils.exceptions import OutputError, FolderCreatorError, FileTypeError
from primer.slice_data import SliceData
from designer.output_data_classes import (
    SlicerOutputData,
    TargetonCSVData,
    ScoringOutputData,
    PrimerDesignerOutputData,
)

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

if TYPE_CHECKING:  # For avoiding circular import dependencies, only import for type checking.
    from src.primer_designer import PrimerDesigner
    from src.cli import Scoring


def timestamped_dir(prefix):
    try:
        FolderCreator.create_timestamped(prefix)
    except FolderCreatorError as err:
        raise OutputError(f'Error creating folder: {err}')
    return FolderCreator.get_dir()


def write_slicer_output(dir_prefix: str, slices: List[dict]) -> SlicerOutputData:
    export_dir = timestamped_dir(dir_prefix)
    result = SlicerOutputData(export_dir)

    result.bed = write_slicer_bed_output(export_dir, slices)
    result.fasta = write_slicer_fasta_output(export_dir, slices)

    print('Slice files saved: ', result.bed, result.fasta)

    return result


def write_slicer_bed_output(export_dir: str, slices: List[dict]) -> str:
    BED_OUTPUT = 'slicer_output.bed'

    bed_path = path.join(export_dir, BED_OUTPUT)
    slices.saveas(bed_path)

    return bed_path


def write_slicer_fasta_output(export_dir: str, slices: List[dict]) -> str:
    FASTA_OUTPUT = 'slicer_output.fasta'

    fasta_path = path.join(export_dir, FASTA_OUTPUT)
    slices.save_seqs(fasta_path)

    return fasta_path


def export_to_bed(bed_rows: list, export_dir: str) -> str:
    PRIMER_OUTPUT_BED = 'p3_output.bed'

    p3_bed = BedTool(bed_rows)
    bed_path = path.join(export_dir, PRIMER_OUTPUT_BED)
    p3_bed.saveas(bed_path)

    return bed_path


def export_to_csv(data: Union[list, dict], export_dir: str, filename: str,
        headers: List[str], delimiter: str = ',') -> str:
    kwargs = {'delimiter': delimiter, 'fieldnames': headers}
    csv_path = Path(export_dir) / filename

    with open(csv_path, "w", newline='') as f:
        output_writer = csv.DictWriter(f, **kwargs)
        output_writer.writeheader()
        if isinstance(data, dict):
            output_writer.writerow(data)
        else:
            output_writer.writerows(data)

    return csv_path


def write_targeton_csv(
        ipcress_input: str,
        slices: List[SliceData],
        dirname: str,
        dir_timestamped=False
) -> TargetonCSVData:
    TARGETON_CSV = 'targetons.csv'

    csv_rows = []
    with open(ipcress_input) as fh:
        ipcress_input_data = fh.read()
    for slice in slices:
        # corresponding primer pair names will be prefixed by region name
        primer_pair_iterator = re.finditer(rf'^{slice.name}\S*', ipcress_input_data, re.MULTILINE)
        for primer_pair in primer_pair_iterator:
            csv_rows.append([primer_pair.group(), slice.name])

    if not dir_timestamped:
        dirname = timestamped_dir(dirname)
    csv_path = path.join(dirname, TARGETON_CSV)
    with open(csv_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerows(csv_rows)

    print(f'Targeton csv generated: {csv_path}')

    result = TargetonCSVData(dirname)
    result.csv = csv_path
    return result


def write_scoring_output(scoring: Scoring, output_tsv: str) -> ScoringOutputData:
    scoring.save_mismatches(output_tsv)

    result = ScoringOutputData('')
    result.tsv = output_tsv

    print(f'Scoring file saved: {output_tsv}')

    return result


def export_primer_design_to_file(primer_designer: PrimerDesigner, filename: str, export_dir: str, file_type: str) -> str:
    accepted_file_types = [r'.json', r'.csv']
    if file_type not in accepted_file_types:
        raise FileTypeError(f"Unknown filetype passed {file_type}.")

    filename = Path(filename)
    if not filename.suffix:
        filename = filename.with_suffix(file_type)
    path = export_dir / filename
    with open(path, 'w') as f:
        if file_type == r'.json':
            primer_designer.dump_json(f, sort_keys=True, indent=4)
        elif file_type == r'.csv':
            flat_dict_list = primer_designer.flatten()
            writer = csv.DictWriter(f, fieldnames=list(flat_dict_list[0].keys()))
            writer.writeheader()
            writer.writerows(flat_dict_list)

    return str(path)


def write_primer_design_output(
    primer_designer : PrimerDesigner,
    prefix='',
    existing_dir='',
) -> PrimerDesignerOutputData:

    if existing_dir:
        export_dir = existing_dir
    else:
        export_dir = timestamped_dir(prefix)

    result = PrimerDesignerOutputData(export_dir)
    filename = r'primer_designer'
    result.csv = export_primer_design_to_file(primer_designer, filename, export_dir, '.csv')
    result.json = export_primer_design_to_file(primer_designer, filename, export_dir, '.json')
    result.dir = export_dir
    print(f'Primer Designer files saved:{result.csv}, {result.json}')

    return result

def export_retrieved_fasta(slice_data, export_dir: str) -> str:

    if not slice_data.bases:
        raise ValueError("Retrieved sequence is empty")


    if slice_data.strand not in {"+", "-"}:
        raise ValueError(f"Invalid strand value: {slice_data.strand}")

    os.makedirs(export_dir, exist_ok=True)

    targeton_id = slice_data.name
    filename = f"{targeton_id}_retrieved.fa"
    fasta_path = path.join(export_dir, filename)

    header_id = (
                    f"{targeton_id}:extended:GRCh38:"
                    f"{slice_data.chromosome}:"
                    f"{slice_data.start}-{slice_data.end}"
                    f"({slice_data.strand}):"
                    f"{slice_data.flanking_region}"
                )
    sequence = slice_data.bases

    record = SeqRecord(
                    Seq(sequence),
                    id=header_id,
                    description="")


    try:
        SeqIO.write(record, fasta_path, "fasta")
    except Exception as e:
        raise IOError(f"Failed to write FASTA file: {fasta_path}") from e

    return fasta_path
