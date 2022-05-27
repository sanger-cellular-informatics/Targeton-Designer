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
    parser.add_argument('bed',
        help='BED file containing regions of interest')
    parser.add_argument('-f', '--fasta',
        help='FASTA file to retrieve sequences from',
        default='')
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
    return parser.parse_args(args)

def main(params):
    print(params)

    bed_path = params['bed']
    bed_file = os.path.basename(bed_path)
    bed_param = f'/tmp/{bed_file}'

    fasta_param = ''
    fasta_path = params['fasta']
    if fasta_path:
        fasta_file = os.path.basename(fasta_path)
        fasta_param = f'--fasta /tmp/{fasta_file}'

    subprocess.run(['cp', bed_path, f'./{bed_file}'])
    
    if fasta_path:
        subprocess.run(['cp', fasta_path, f'./{fasta_file}'])
    
    pwd = subprocess.run(['pwd'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
    cmd = f'docker run -i -t --rm -v targeton-designer_data:/data -v {pwd}/:/tmp/ targeton-designer_slicer {bed_param} {fasta_param}'
    
    os.system(cmd)
    
    return 1

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    print(main(vars(args)))
