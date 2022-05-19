import primer3
import json
import re
import os
import collections

def primer3_runner(design_input, strand):
    print('Designing primers for the region')
    design = primer3_design(design_input)
    print('Pick primers')
    primers = locate_primers(design, design_input['SEQUENCE_ID'], strand)
    print(primers)
    return primers

def primer3_design(primer3_input):
    p3_config_loc = os.environ.get('PRIMER3_CONFIG')
    
    primer3_config = {}
    with open(os.path.join(os.path.dirname(__file__), p3_config_loc), "r") as p3:
        primer3_config = json.load(p3)

    design = primer3.bindings.designPrimers(primer3_input, primer3_config)

    return design

def locate_primers(design, slice_name, strand):
    primer_keys = design.keys()
    primers = collections.defaultdict(dict)
    for key in primer_keys:
        primer_details = capture_primer_details(key)
        if primer_details:
            primer_id = primer_details['id']
            primer_field = primer_details['field']
            pair_number = primer_details['pair']
            
            libamp_name = name_primers(primer_details, strand)
            
            primer_name = slice_name + "_" + libamp_name + "_" + pair_number
            primers[primer_name][primer_field] = design[key]
    
    return primers

def capture_primer_details(primer_name):
    match = re.search(r'^(primer_(left|right)_(\d+))\_(\S+)$', primer_name.lower())
    result = {}
    if match:
        primer_id = match.group(1)
        primer_side = match.group(2)
        pair_number = match.group(3)
        primer_field = match.group(4)

        result = {
            'id' : primer_id,
            'side' : primer_side,
            'field' : primer_field,
            'pair' : pair_number
        }

    return result

def name_primers(primer_details, strand):
    fwd_primers = {
        'left' : 'LibAmpF',
        'right' : 'LibAmpR',
    }
    rev_primers = {
        'left' : 'LibAmp_R',
        'right' : 'LibAmp_F',
    }
    names = {
        '1' : fwd_primers,
        '-1' : rev_primers,
    }

    primer_name = names[strand][primer_details['side']]

    return primer_name

