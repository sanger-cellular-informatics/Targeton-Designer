#!/usr/bin/env python3
from pybedtools import BedTool
import argparse
import sys

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

def positive_int(arg):
    if int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0')
    return int(arg)

def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Get sequence slices for regions in'\
            ' BED file according to parameters specified')
    parser.add_argument('input_bed',
        help='BED file containing regions of interest')
    parser.add_argument('input_fasta',
        help='FASTA file to retrieve sequences from')
    parser.add_argument('-f5', '--flank_5',
        help='how far to extend region at 5\' end (default 50nt)',
        type=int, default=50)
    parser.add_argument('-f3', '--flank_3',
        help='how far to extend region at 3\' end (default 50nt)',
        type=int, default=50)
    parser.add_argument('-l', '--length',
        help='length of each slice (default 210nt)',
        type=positive_int, default=210)
    parser.add_argument('-o', '--offset',
        help='offset between each slice (default 5nt)',
        type=positive_int, default=5)
    parser.add_argument('--output_fasta',
        help='output slice sequences to fasta file')
    parser.add_argument('--output_bed',
        help='output bed file with slice coordinates')
    return parser.parse_args(args)

def get_slices(params):
    bed = BedTool(params['input_bed'])
    slice_bed = BedTool(get_slice_data(bed, params))
    # return named slice sequences on specified strand
    return slice_bed.sequence(fi=params['input_fasta'], name=True, s=True)

def main(params):
    slices = get_slices(params)
    if params['output_bed'] is not None:
        slices.saveas(params['output_bed'])
    if params['output_fasta'] is not None:
        slices.save_seqs(params['output_fasta'])
        print('Slice sequences saved!')
    else:
        print(slices.print_sequence())

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    main(vars(args))
