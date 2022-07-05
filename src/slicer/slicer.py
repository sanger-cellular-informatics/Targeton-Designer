#!/usr/bin/env python3
import re
from pybedtools import BedTool
from pybedtools.helpers import BEDToolsError

class SlicerError(Exception):
    pass

def handle_one_based_input(input_bed):
    adjusted_tsv = []
    with open(input_bed) as file:
        tsv = csv.reader(file, delimiter="\t")
        adjusted_tsv = decrement_one_based_starts(tsv, adjusted_tsv)
    return adjusted_tsv

def decrement_one_based_starts(tsv, new_tsv):
    #BED is only 0-based on the start thus only need to edit column 1
    for row in tsv:
        row[1] = str(int(row[1]) - 1)
        new_tsv.append(row)
    return new_tsv

def _generate_slice_data(exon, exon_name, params):
    slices = []
    start = exon.start - params['flank_5']
    end = start + params['length']
    count = 1
    while end <= (exon.end + params['flank_3']):
        slice_name = f'{exon_name}_{count}'
        slices.append((exon.chrom, start, end, slice_name,
            exon.score, exon.strand))
        start += params['offset']
        end += params['offset']
        count += 1
    return slices

def get_slice_data(bed, params):
    slices = []
    count = 1
    for exon in bed:
        name = exon.name if exon.name != '.' else f'region{count}'
        slices.extend(_generate_slice_data(exon, name, params))
        count += 1
    return slices

def get_slices(params):

    input_bed = params['bed']
    if params['1b']:
        input_bed = handle_one_based_input(params['bed'])
    bed = BedTool(input_bed)

    slice_bed = BedTool(get_slice_data(bed, params))
    # return named, coords slice sequences on specified strand
    seq_options = {
        "fi"    : params['fasta'],
        "s"     : True,
        "name+" : True
    }
    seq = {}
    try:
        seq = slice_bed.sequence(**seq_options)
    except BEDToolsError as bed_err:
        if not re.search(r'\*{5}ERROR:\ Unrecognized parameter: -name\+\ \*{5}', bed_err.args[1]):
            template = "PyBEDTools exited with err type {0}. Arguments:\n{1!r}"
            message = template.format(type(bed_err).__name__, bed_err.args[1])
            raise BEDToolsError(bed_err, message)
        del seq_options['name+']
        seq_options['name'] = True
        seq = slice_bed.sequence(**seq_options)

    return seq

def main(params):
    try:
        slices = get_slices(params)
        return slices

    except Exception as err:
        raise SlicerError('Unexpected error occurred: {0}'.format(err))

if __name__ == '__main__':
    main()