import os
import argparse
import sys
import subprocess


def positive_int(arg):
    if int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0')
    return int(arg)

def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Get sequence slices for regions in'\
            ' BED file according to parameters specified')
    parser.add_argument('-b', '--bed',
        help='BED file containing regions of interest')
    parser.add_argument('-f', '--fasta',
        help='FASTA file to retrieve sequences from')
    parser.add_argument('-f5', '--flank_5',
        help='how far to extend region at 5\' end',
        type=int, default=50)
    parser.add_argument('-f3', '--flank_3',
        help='how far to extend region at 3\' end',
        type=int, default=50)
    parser.add_argument('-l', '--length',
        help='length of each slice',
        type=positive_int, default=210)
    parser.add_argument('-o', '--offset',
        help='offset between each slice',
        type=positive_int, default=5)
    return parser.parse_args(args)

def main(params):

    bed_path = params['bed']
    bed_file = os.path.basename(bed_path)
    bed_param = f'/tmp/{bed_file}'

    subprocess.run(['mkdir', '/tmp/slicer'])

    subprocess.run(['cp', bed_path, f'/tmp/slicer/{bed_file}'])

    fasta_path = params['fasta']
    if fasta_path:
        fasta_file = os.path.basename(fasta_path)
        fasta_param = f'/tmp/{fasta_file}'
        subprocess.run(['cp', fasta_path, f'/tmp/slicer/{fasta_file}'])
    else:
        fasta_param = '/data/reference.fa'

    cmd = ['docker', 'run', '-i', '-t', '--rm',
           '-v', 'targetondesigner_data:/data',
           '-v', '/tmp/slicer:/tmp/',
           'targetondesigner_slicer',
           '-b', bed_param, 
           '-f', fasta_param,
           '--flank_5', str(params['flank_5']),
           '--flank_3', str(params['flank_3']),
           '--length', str(params['length']),
           '--offset', str(params['offset'])]
    print(cmd)
    output = subprocess.run(cmd,
        stdout=subprocess.PIPE).stdout.decode('utf-8').strip()

    subprocess.run(['rm', '-r', '/tmp/slicer'])
    
    return output

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    print(main(vars(args)))
