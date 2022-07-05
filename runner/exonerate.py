import sys
import collections
import re
import os
import subprocess
import argparse
import csv

from os import path
from pprint import pprint

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
    return parser.parse_args(args)

def run_ipcress(run_id, dir_path, reference_file):
    primers = retrieve_p3_output(dir_path)
    pprint(primers)
    input_path = generate_ipcress_input(run_id, dir_path)

    input_cmd = '--input ' + input_path
    seq = ' --sequence ' + reference_file
    mismatch = ' --mismatch 5'
    cmd = "ipcress " + input_cmd + seq + mismatch
    ipcress = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    
    stnd, err = ipcress.communicate()
    print("stdout:", stnd)
    print("stderr:", err)

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
        # ENSE00003571441_HG6_6_LibAmpR_0
        match = re.search(r'^(\w+_LibAmp)(F|R)_(\d+)$', row['primer'])
        if match:
            key = match.group(1) + '_' + match.group(3)
            primer_data[key][match.group(2)] = row['sequence']

    return primer_data







def generate_ipcress_input(run_id, dir_path):
    primers = {}
    pairs = pair_primers(primers)
    rows = format_ipcress_input(run_id, pairs)
    file_path = write_ipcress_input(run_id, dir_path, rows)

    return file_path    

def write_ipcress_input(run_id, dir_path, rows):
    path = dir_path + '/' + run_id + '/' + 'primer_input.txt'
   
    os.makedirs(os.path.dirname(path), exist_ok=True)    
    file_h = open(path, "w")
    for row in rows:
        file_h.write(row + "\n")
    file_h.close

    return path

def format_ipcress_input(run_id, pairs):
    ipcress_input = []
    rows = pairs.keys()
    
    for key in rows:
        row_id = run_id + '_' + key
        left = pairs[key]['left']['sequence']
        right = pairs[key]['right']['sequence']
        min_val = '200'
        max_val = '400'
        line = ' '.join([
            row_id,
            left,
            right,
            min_val,
            max_val
        ])
        ipcress_input.append(line)
    return ipcress_input    

def pair_primers(primers):
    primer_keys = primers.keys()
    primers_prep = collections.defaultdict(dict)
    for key in primer_keys:
        match = re.search(r'^primer_(left|right)_(\d+)$', key)
        if match:
            primer_side = match.group(1)
            primer_pair = match.group(2)
            primers_prep[primer_pair][primer_side] = primers[key]
    return primers_prep

def main(params):
    run_ipcress('test', params['dir'], params['ref'])

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    main(vars(args))
