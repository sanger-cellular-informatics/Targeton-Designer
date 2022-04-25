import docker
import collections
import re
import os

def run_exonerate_docker():
    client = docker.from_env()
    test = client.containers.run('exonerate','--version')
    print(test, flush=True)

def run_ipcress(run_id, primers, dir_path):
    input_path = generate_ipcress_input(run_id, primers, dir_path)

def generate_ipcress_input(run_id, primers, dir_path):
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
