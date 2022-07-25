import sys
from collections import defaultdict
import re
import os
import subprocess
import argparse
import csv

from os import path
from dataclasses import dataclass

from utils.file_system import write_to_text_file, read_csv_to_dict

@dataclass
class IpcressResult:
    stnd: bytes
    err: bytes

def ipcress_runner(params) -> IpcressResult:
    if params['primers']:
        print('Loading custom iPCRess input file')
        input_path = params['primers']
        result = run_ipcress(input_path, params)
    else:
        print('Building iPCRess input file.')
        primer_data = extract_primer_data(params['p3_csv'])
        input_path = get_ipcress_input(primer_data, params)
        result = run_ipcress(input_path, params)
        validate_primers(result.stnd.decode(),
            primer_data, params['pretty'])

    print('Finished!')

    return result

def run_ipcress(input_path, params) -> IpcressResult:
    cmd = ' '.join(['ipcress', input_path, params['fasta'],
        '--mismatch', params['mismatch']])

    cmd = prettify_output(params['pretty'], cmd)

    print('Running Exonerate iPCRess with the following command:')
    print(cmd)
    
    ipcress = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    stnd, err = ipcress.communicate()
    return IpcressResult(stnd, err)

def prettify_output(prettify, cmd):
    pretty_opt = 'false'
    if prettify:
        pretty_opt = 'true'

    return cmd + ' --pretty ' + pretty_opt

def extract_primer_data(p3_csv):
    csv_obj = read_csv_to_dict(p3_csv)
    primer_data = defaultdict(lambda: defaultdict(dict))
    for row in csv_obj:
        #Capture primer name and orientation
        #ENSE00003571441_HG6_6_LibAmpR_0
        match = re.search(r'^(\w+_LibAmp)(F|R)_(\d+)$', row['primer'])
        if match:
            key = match.group(1) + '_' + match.group(3)
            primer_data[key][match.group(2)]['seq'] = row['sequence']
            primer_data[key][match.group(2)]['start'] = row['primer_start']

    return primer_data

def get_ipcress_input(primer_data, params):
    formatted_primers = format_ipcress_primers(params['min'],
        params['max'], primer_data)

    input_path = write_to_text_file(params['dir'],
        formatted_primers, 'ipcress_primer_input')

    return input_path

def format_ipcress_primers(min_amp, max_amp, primers):
    ipcress_input = []
    rows = primers.keys()
    
    for key in rows:
        left = primers[key]['F']['seq']
        right = primers[key]['R']['seq']

        line = ' '.join([
            key,
            left,
            right,
            min_amp,
            max_amp
        ])
        ipcress_input.append(line)

    return ipcress_input    

def validate_primers(ipcress_output, primer_data, pretty):
    if pretty:
        print('Output is pretty, skipping validation')
        return
    print('Validating primers...')
    for primer_pair in primer_data.keys():
        fwd_coord = primer_data[primer_pair]['F']['start']
        rev_coord = primer_data[primer_pair]['R']['start']
        match = re.search((
            fr'{primer_pair} \d+ A {fwd_coord} 0 '
            fr'B {rev_coord} 0 forward$'), ipcress_output)
        if not match:
            print(f'No valid primer pair found for {primer_pair}')

def main(params):
    return ipcress_runner(params)

if __name__ == '__main__':
    main()
