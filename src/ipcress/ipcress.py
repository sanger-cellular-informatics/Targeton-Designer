import sys
import collections
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
    stnd: str
    err: str

def run_ipcress(params) -> IpcressResult:
    input_path = determine_ipcress_input(params, primers_txt=params['primers'])
    cmd = "ipcress " + input_path + ' ' + params['fasta'] + ' --mismatch ' + params['mismatch']

    cmd = prettify_output(params['pretty'], cmd)

    print('Running Exonerate iPCRess with the following command:')
    print(cmd)
    
    ipcress = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    stnd, err = ipcress.communicate()
    result = IpcressResult(stnd, err)
    
    print('Finished!')

    return result

def prettify_output(prettify, cmd):
    pretty_opt = 'false'
    if prettify == 'true':
        pretty_opt = 'true'

    return cmd + ' --pretty ' + pretty_opt

def determine_ipcress_input(params, primers_txt = ''):
    input_path = ''
    if primers_txt:
        print('Loading custom iPCRess input file')
        input_path = primers_txt
    else:
        print('Building iPCRess input file.')
        input_path = retrieve_primer3_output(params)
    return input_path

def format_ipcress_primers(min_amp, max_amp, primers):
    ipcress_input = []
    rows = primers.keys()
    
    for key in rows:
        left = primers[key]['F']
        right = primers[key]['R']

        line = ' '.join([
            key,
            left,
            right,
            min_amp,
            max_amp
        ])
        ipcress_input.append(line)

    return ipcress_input    

def retrieve_primer3_output(params):
    file_data = read_csv_to_dict(params['p3_csv'])
    
    primer_data = extract_primer_sequences(file_data)
    formatted_primers = format_ipcress_primers(params['min'], params['max'], primer_data)
    
    input_path = write_to_text_file(params['dir'], formatted_primers, 'ipcress_primer_input')
    
    return input_path
    
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
    return run_ipcress(params)

if __name__ == '__main__':
    main()
