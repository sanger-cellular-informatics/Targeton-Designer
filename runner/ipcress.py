import sys
import collections
import re
import os
import subprocess
import argparse
import csv

from os import path

def add_modules_to_sys_path():
    BASE_PATH = path.dirname(path.dirname(path.abspath(__file__)))
    sys.path.append(BASE_PATH)

add_modules_to_sys_path()

from utils.file_system import check_file_exists

def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Primer scoring using Exonerate iPCRess')
    parser.add_argument('--ref',
        help='Genomic Reference File.')
    parser.add_argument('--dir',
        help='Output directory created by slicer tool')
    parser.add_argument('--min',
        help='Minimum amplicon length',
        default='200')
    parser.add_argument('--max',
        help='Maximum amplicon length',
        default='300')
    return parser.parse_args(args)

def run_ipcress(run_id, params):
    dir_path = params['dir']
    reference_file = params['ref']
    
    print('Building iPCRess input file.')

    primers = retrieve_p3_output(dir_path)
    formatted_primers = format_ipcress_primers(primers, params)
    input_path = write_ipcress_input_file(dir_path, formatted_primers)

    input_cmd = '--input ' + input_path
    seq = ' --sequence ' + reference_file
    mismatch = ' --mismatch 5'
    cmd = "ipcress " + input_cmd + seq + mismatch

    print('Running Exonerate iPCRess.')
    
    ipcress = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    
    stnd, err = ipcress.communicate()

    print('Finished!')
    print("stdout:", stnd)
    print("stderr:", err)

def format_ipcress_primers(primers, params):
    ipcress_input = []
    rows = primers.keys()
    
    for key in rows:
        left = primers[key]['F']
        right = primers[key]['R']

        line = ' '.join([
            key,
            left,
            right,
            params['min'],
            params['max']
        ])
        ipcress_input.append(line)

    return ipcress_input    

def write_ipcress_input_file(dir_path, rows):
    path = dir_path + '/' + 'ipcress_primer_input.txt'
   
    file_h = open(path, "w")
    for row in rows:
        file_h.write(row + "\n")
    file_h.close

    return path

def retrieve_p3_output(dir_path):
    p3_csv = path.join(dir_path, 'p3_output.csv')
    check_file_exists(p3_csv)
    
    primer_data = {}
    with open(p3_csv) as csv_file:
        csv_obj = csv.DictReader(csv_file, delimiter=',')
        primer_data = extract_primer_sequences(csv_obj)
    
    return primer_data
    
def extract_primer_sequences(csv_obj):
    primer_data = collections.defaultdict(dict)
    for row in csv_obj:
        #Capture primer name and orientation
        #ENSE00003571441_HG6_6_LibAmpR_0
        match = re.search(r'^(\w+_LibAmp)(F|R)_(\d+)$', row['primer'])
        if match:
            key = match.group(1) + '_' + match.group(3)
            primer_data[key][match.group(2)] = row['sequence']

    return primer_data

def main(params):
    run_ipcress('test', params)

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    main(vars(args))
