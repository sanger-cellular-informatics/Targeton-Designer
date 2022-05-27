import os
import argparse
import sys
import subprocess


def positive_int(arg):
    if int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0')
    return int(arg)


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('--seqs', action='store', type=str, nargs='+', help='Path to sequences file')
    parser.add_argument('--dir', action='store', type=str, nargs='+', help='Output folder location')
    parser.add_argument('--ref', action='store', type=str, nargs='+', help='Reference file path')
    return parser.parse_args(args)


def main(params):
    print(params)

    seq_path = params.seqs[0]
    seq_file = os.path.basename(seq_path)
    seqs_param = f'--seqs /tmp/{seq_file}'

    dir_path = params.dir[0]
    dir_file = os.path.basename(dir_path)
    dir_param = f'--dir /tmp/{dir_file}'

    ref_path = params.ref[0]
    ref_file = os.path.basename(ref_path)
    ref_param = f'--ref /tmp/{ref_file}'

    subprocess.run(['mkdir', '/tmp/primer3/'])

    subprocess.run(['cp', seq_path, f'/tmp/primer3/{seq_file}'])
    subprocess.run(['cp', dir_path, f'/tmp/primer3/{dir_file}'])
    subprocess.run(['cp', ref_path, f'/tmp/primer3{ref_file}'])

    #pwd = subprocess.run(['pwd'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
    cmd = f'docker run -i -t --rm -v targeton-designer_data:/data -v /tmp/primer3/:/tmp/ targeton-designer_primer3 {seqs_param} {dir_param} {ref_param}'

    os.system(cmd)

    subprocess.run(['rm', '-r', '/tmp/primer3'])

    return 1


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    print(main(args))
