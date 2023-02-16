import csv
from Bio import SeqIO
import re

from utils.exceptions import FileFormatError, FileValidationError
from utils.file_system import check_file_exists


def validate_bed_format(bed:str):
    with open(bed, newline='') as file:
        tsv_file = csv.reader(file, delimiter='\t')

        line_num = 1
        for line in tsv_file:
            if len(line) < 6:
                raise FileFormatError(f'Unable to read in BED file correctly. Check file format on line {line_num}.')

            line_num = line_num + 1


def validate_fasta_format(fasta:str):
    with open(fasta) as handle:
        if not any(SeqIO.parse(handle, "fasta")):
            raise FileFormatError('Unable to read in FastA file correctly. Check file format.')


def validate_bed_content(bed:str):
    with open(bed, newline='') as file:
        tsv_file = csv.reader(file, delimiter='\t')
        line_num = 1
        for line in tsv_file:

            if not re.search(r'^(?:[Cc][Hh][Rr])?([XxYy]|[1-9]|1\d|2[012])$', line[0]):
                raise ValueError(f'Chromosome format incorrect on line {line_num}: {line[0]}')

            if not re.search(r'^\d+$', line[1]):
                raise ValueError(f'Start coordinate format incorrect on line {line_num}: {line[1]}')

            if not re.search(r'^\d+$', line[2]):
                raise ValueError(f'End coordinate format incorrect on line {line_num}: {line[2]}')

            if int(line[2]) < int(line[1]):
                raise ValueError(f'End coordinate must be greater than start coordinate '
                                 f'on line {line_num}. Start: {line[1]} End: {line[2]}')

            if (int(line[2]) - int(line[1])) > 10000:
                raise ValueError(f'Difference between start coordinate and end coordinate '
                                 f'must be less than 10000. On line {line_num} '
                                 f'Difference: {int(line[2]) - int(line[1])}')

            if not line[3]:
                raise ValueError(f'Error with name field, if no name is supplied please mark '
                                 f'with a \'.\' on line {line_num}: {line[3]}')

            if not line[4]:
                raise ValueError(f'Error with score field, if no score is supplied please mark '
                                 f'with a \'.\' on line {line_num}: {line[4]}')

            if not re.search(r'^[+-]$', line[5]):
                raise ValueError(f'Strand format incorrect on line {line_num}: {line[5]}')

            line_num = line_num + 1


def validate_p3_csv(p3_csv:str):
    with open(p3_csv, newline='') as csv_file:
        data = csv.DictReader(csv_file, delimiter=',')
        expected_cols = [
            'primer', 'sequence', 'chr', 'primer_start', 'primer_end',
            'tm', 'gc_percent', 'penalty', 'self_any_th', 'self_end_th',
            'hairpin_th', 'end_stability'
        ]
        if check_if_missing_fields(data, expected_cols):
            raise FileFormatError(f'Missing columns in Primer3 CSV')


def validate_score_tsv(tsv:str):
    with open(tsv, newline='') as tsv_file:
        data = csv.DictReader(tsv_file, delimiter='\t')
        expected_cols = ['Targeton', 'Primer pair', 'A/B/Total', 'WGE format', 'Score']
        if check_if_missing_fields(data, expected_cols):
            raise FileFormatError(f'Missing columns in Scoring TSV')


def validate_files(bed='', fasta='', txt='', p3_csv='', score_tsv=''):
    try:
        if bed:
            check_file_exists(bed)
            validate_bed_format(bed)
            validate_bed_content(bed)

        if fasta:
            check_file_exists(fasta)
            validate_fasta_format(fasta)

        if txt:
            check_file_exists(txt)

        if p3_csv:
            check_file_exists(p3_csv)
            validate_p3_csv(p3_csv)

        if score_tsv:
            check_file_exists(score_tsv)
            validate_score_tsv(score_tsv)

    except ValueError as valErr:
        print('Error occurred while checking file content: {0}'.format(valErr))
    except FileFormatError as fileErr:
        print('Error occurred while checking file format: {0}'.format(fileErr))
    except FileNotFoundError as fileErr:
        print('Input file not found: {0}'.format(fileErr))
    except Exception as err:
        print('Unexpected error occurred: {0}'.format(err))

    return


def check_if_missing_fields(data: dict, fields: list) -> bool:
    missing_fields = []
    for field in fields:
        if field not in data.fieldnames:
            missing_fields.append(field)
    if any(missing_fields):
        return True
    return False
