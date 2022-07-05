import csv
from Bio import SeqIO
import re

from .exceptions import FileFormatError, FileValidationError
from .file_system import check_file_exists

def validate_bed_format(bed):
    with open(bed) as file:
        tsv_file = csv.reader(file, delimiter='\t')

        line_num = 1
        for line in tsv_file:
            if len(line) < 6:
                raise FileFormatError(f'Unable to read in BED file correctly. Check file format on line {line_num}.')

            line_num = line_num + 1

def validate_fasta_format(fasta):
    with open(fasta) as handle:
        if not any(SeqIO.parse(handle, "fasta")):
            raise FileFormatError('Unable to read in Fasta file correctly. Check file format.')

def validate_bed_content(bed):
    with open(bed) as file:
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
                raise ValueError(f'Error with  name field, if no name is supplied please mark '
                                 f'with a \'.\' on line {line_num}: {line[3]}')

            if not line[4]:
                raise ValueError(f'Error with score field, if no score is supplied please mark '
                                 f'with a \'.\' on line {line_num}: {line[4]}')

            if not re.search(r'^[+-]$', line[5]):
                raise ValueError(f'Strand format incorrect on line {line_num}: {line[5]}')

            line_num = line_num + 1

def validate_files(bed, fasta):
    try:
        check_file_exists(bed)
        check_file_exists(fasta)

        validate_bed_format(bed)
        validate_fasta_format(fasta)

        validate_bed_content(bed)

    except ValueError as valErr:
        raise SlicerError('Error occurred while checking file content: {0}'.format(valErr))
    except FileFormatError as fileErr:
        raise SlicerError('Error occurred while checking file format: {0}'.format(fileErr))
    except FileNotFoundError as fileErr:
        raise FileValidationError('Input file not found: {0}'.format(fileErr))
    except Exception as err:
        raise FileValidationError('Unexpected error occurred: {0}'.format(err))

    return