from __future__ import annotations

import csv
import re

from typing import TYPE_CHECKING, Any, List, Tuple
from os import path
from pathlib import Path
from dataclasses import dataclass
from pybedtools import BedTool
from utils.file_system import write_to_text_file, FolderCreator
from utils.exceptions import OutputError, FolderCreatorError, FileTypeError
if TYPE_CHECKING: # For avoiding circular import dependencies, only import for type checking.
    from src.primer_designer import PrimerDesigner
    from src.cli import Scoring


@dataclass
class OutputFilesData:
    dir: str

    def fields(self):
        return list(self.__dataclass_fields__.keys())


@dataclass
class SlicerOutputData(OutputFilesData):
    bed: str = ''
    fasta: str = ''


@dataclass
class PrimerOutputData(OutputFilesData):
    bed: str = ''
    csv: str = ''


@dataclass
class PrimerDesignerOutputData(OutputFilesData):
    csv: str = ''
    json: str = ''


@dataclass
class IpcressOutputData(OutputFilesData):
    input_file: str = ''
    stnd: str = ''
    err: str = ''


@dataclass
class TargetonCSVData(OutputFilesData):
    csv: str = ''


@dataclass
class ScoringOutputData(OutputFilesData):
    tsv: str = ''


@dataclass
class DesignOutputData(OutputFilesData):
    slice_bed: str = ''
    slice_fasta: str = ''
    p3_bed: str = ''
    p3_csv: str = ''
    pd_csv: str = ''
    pd_json: str = ''
    ipcress_input: str = ''
    ipcress_output: str = ''
    ipcress_err: str = ''
    targeton_csv: str = ''
    scoring_tsv: str = ''


def timestamped_dir(prefix):
    try:
        FolderCreator.create_timestamped(prefix)
    except FolderCreatorError as err:
        raise OutputError(f'Error creating folder: {err}')
    return FolderCreator.get_dir()


def write_slicer_output(dir_prefix:str, slices:List[dict]) -> SlicerOutputData:
    export_dir = timestamped_dir(dir_prefix)
    result = SlicerOutputData(export_dir)

    result.bed = write_slicer_bed_output(export_dir, slices)
    result.fasta = write_slicer_fasta_output(export_dir, slices)

    print('Slice files saved: ', result.bed, result.fasta)

    return result


def write_slicer_bed_output(export_dir:str, slices:List[dict]) -> str:
    BED_OUTPUT = 'slicer_output.bed'

    bed_path = path.join(export_dir, BED_OUTPUT)
    slices.saveas(bed_path)

    return bed_path


def write_slicer_fasta_output(export_dir:str, slices:List[dict]) -> str:
    FASTA_OUTPUT = 'slicer_output.fasta'

    fasta_path = path.join(export_dir, FASTA_OUTPUT)
    slices.save_seqs(fasta_path)

    return fasta_path


def export_primers_to_csv(slices:List[dict], export_dir:str) -> str:
    PRIMER3_OUTPUT_CSV = 'p3_output.csv'

    headers = ['primer', 'sequence', 'chr', 'primer_start', 'primer_end', 'tm', 'gc_percent',
               'penalty', 'self_any_th', 'self_end_th', 'hairpin_th', 'end_stability']
    rows = construct_csv_format(slices, headers)

    csv_path = export_to_csv(rows, export_dir, PRIMER3_OUTPUT_CSV, headers, delimiter=',')
    return csv_path

def export_to_csv(data:Tuple[list, dict], export_dir:str, filename:str, headers:List[str], delimiter=',') -> str:
    writers = {dict:csv.DictWriter, list:csv.writer}
    kwargs = {'delimiter':delimiter}
    csv_path = Path(export_dir)/filename
    if isinstance(data, dict):
        kwargs['fieldnames'] = headers
        writer = writers[dict]
    else:
        writer = writers[list]

    with open(csv_path, "w", newline='') as f:
        output_writer = writer(f, **kwargs)
        if isinstance(data, dict):
            output_writer.writeheader()
        output_writer.writerows(data)

    return csv_path

def construct_csv_format(slices:List[dict], headers:list) -> list:
    rows = []

    for slice_data in slices:
        primers = slice_data['primers']
        for primer in primers:
            primers[primer]['primer'] = primer
            primers[primer]['chr'] = slice_data['chrom']

            del primers[primer]['coords']
            del primers[primer]['side']
            del primers[primer]['strand']

            rows.append(primers[primer])

    return rows


def construct_bed_format(slices:List[dict])-> list:
    rows = []
    for slice_data in slices:
        primers = slice_data['primers']
        for primer in primers:
            primer_data = primers[primer]
            # chr,chrStart,chrEnd,name,score,strand
            # Score unknown until iPCRess
            row = [
                slice_data['chrom'],
                primer_data['primer_start'],
                primer_data['primer_end'],
                primer,
                '0',
                primer_data['strand']
            ]
            rows.append(row)
    return rows


def export_to_bed(bed_rows:list, export_dir:str) -> str:
    PRIMER_OUTPUT_BED = 'p3_output.bed'

    p3_bed = BedTool(bed_rows)
    bed_path = path.join(export_dir, PRIMER_OUTPUT_BED)
    p3_bed.saveas(bed_path)

    return bed_path


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


def write_ipcress_input(export_dir, formatted_primers) -> str:
    INPUT_FILE_NAME = 'ipcress_primer_input'

    file_path = write_to_text_file(export_dir, formatted_primers, INPUT_FILE_NAME)

    return file_path


def write_ipcress_output(stnd='', err='', existing_dir='') -> IpcressOutputData:
    IPCRESS_OUTPUT_TXT = 'ipcress_output'

    result = IpcressOutputData(existing_dir)

    result.stnd = write_to_text_file(existing_dir, stnd, IPCRESS_OUTPUT_TXT)
    result.err = write_to_text_file(existing_dir, err, IPCRESS_OUTPUT_TXT + "_err")

    return result


def write_targeton_csv(ipcress_input:str, bed:str, dirname:str, dir_timestamped=False) -> TargetonCSVData:
    TARGETON_CSV = 'targetons.csv'
    bed = BedTool(bed)
    csv_rows = []
    with open(ipcress_input) as fh:
        ipcress_input_data = fh.read()
    for region in bed:
        # corresponding primer pair names will be prefixed by region name
        primer_pair_iterator = re.finditer(rf'^{region.name}\S*', ipcress_input_data, re.MULTILINE)
        for primer_pair in primer_pair_iterator:
            csv_rows.append([primer_pair.group(), region.name])

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


def write_scoring_output(scoring:Scoring, output_tsv:str) -> ScoringOutputData:
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
