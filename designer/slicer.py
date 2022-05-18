#!/usr/bin/env python3
from pybedtools import BedTool
import argparse
import sys

def _generate_slice_coordinates(exon, count, params):
    slices = []
    name = exon.name if exon.name != '.' else count
    start = exon.start - params['flank_5']
    end = start + params['length']
    while end <= (exon.end + params['flank_3']):
        slices.append((exon.chrom, start, end, name))
        start += params['offset']
        end += params['offset']
    return slices

def get_slice_coordinates(bed, params):
    slices = []
    count = 1
    for exon in bed:
        slices.extend(_generate_slice_coordinates(exon, count, params))
        count += 1
    return slices

def get_slice_sequences(bed, fasta):
    seqs = {}
    results = bed.sequence(fi=fasta, tab=True,
        name=True).print_sequence().strip()
    for row in results.split('\n'):
        name, sequence = row.split('\t')
        seqs[name] = sequence
    return seqs

def positive_int(arg):
    if int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0')
    return int(arg)

def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Get sequence slices for regions in'\
            ' BED file according to parameters specified')
    parser.add_argument('bed',
        help='BED file containing regions of interest')
    parser.add_argument('fasta',
        help='FASTA file to retrieve sequences from')
    parser.add_argument('-f5', '--flank_5',
        help='how far to extend region at 5\' end',
        type=positive_int, default=0)
    parser.add_argument('-f3', '--flank_3',
        help='how far to extend region at 3\' end',
        type=positive_int, default=0)
    parser.add_argument('-l', '--length',
        help='length of each slice',
        type=positive_int, default=210)
    parser.add_argument('-o', '--offset',
        help='offset between each slice',
        type=positive_int, default=5)
    parser.add_argument('--output_slice_bed',
        help='output bed file with slice coordinates',
        nargs='?', const='slices.bed')
    return parser.parse_args(args)

def main(params):
    bed = BedTool(params['bed'])
    slice_bed = BedTool(get_slice_coordinates(bed, params))
    if 'output_slice_bed' in params:
        slice_bed.saveas(params['output_slice_bed'])
    return get_slice_sequences(slice_bed, params['fasta'])

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    print(main(vars(args)))
