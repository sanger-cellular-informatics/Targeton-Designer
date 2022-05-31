import primer3
import json
import re
import os
import collections
import csv

from pybedtools import BedTool

def primer3_runner(design_input, strand):
    slice_start = 1
    print('Designing primers for the region')
    design = primer3_design(design_input)
    print('Naming primers')
    primers = locate_primers(design, design_input['SEQUENCE_ID'], strand, slice_start)
    print('Exporting to BED and CSV')
    export_primers(primers, 'chr1', slice_start, strand)
 
    return primers

def primer3_design(primer3_input):
    p3_config_loc = os.environ.get('PRIMER3_CONFIG')
    
    primer3_config = {}
    with open(os.path.join(os.path.dirname(__file__), p3_config_loc), "r") as p3:
        primer3_config = json.load(p3)

    design = primer3.bindings.designPrimers(primer3_input, primer3_config)

    return design

def locate_primers(design, slice_name, strand, slice_start):
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
            
            if primer_field == 'coords':
                primer_coords = calculate_primer_coords(primer_details['side'], design[key], slice_start)
                primers[primer_name]['primer_start'] = primer_coords[0] 
                primers[primer_name]['primer_end'] = primer_coords[1] 
    
    return primers

def calculate_primer_coords(side, coords, slice_start):
    slice_start = int(slice_start)
    left_flank = {
        'start' : slice_start,
        'end' : slice_start + int(coords[1])
    }

    slice_end = slice_start + int(coords[0])
    right_flank  = {
        'start' : slice_end - int(coords[1]),
        'end' : slice_end,
    }

    slice_coords = {
        'left' : left_flank,
        'right' : right_flank
    }

    start = slice_coords[side]['start']
    end = slice_coords[side]['end']

    return start, end

def capture_primer_details(primer_name):
    match = re.search(r'^(primer_(left|right)_(\d+))(\_(\S+))?$', primer_name.lower())
    result = {}
    if match:
        primer_id = match.group(1)
        primer_side = match.group(2)
        pair_number = match.group(3)
        primer_field = match.group(5)
        if primer_field == None:
            primer_field = 'coords'
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
        'left' : 'LibAmpR',
        'right' : 'LibAmpF',
    }
    names = {
        '1' : fwd_primers,
        '-1' : rev_primers,
    }

    primer_name = names[strand][primer_details['side']]

    return primer_name

def export_primers(primers, chromo, slice_coords, strand):
    bed_rows = construct_bed_format(primers, chromo, slice_coords, strand)
    export_to_bed(bed_rows)
    export_to_csv(primers)

def export_to_csv(primers):
    headers = ['primer', 'sequence', 'tm', 'gc_percent', 'penalty', 'self_any_th', 'self_end_th', 'hairpin_th', 'end_stability']
    rows = construct_csv_format(primers, headers)

    with open('p3_output.csv', "w") as p3_fh:
        p3_out = csv.DictWriter(p3_fh, fieldnames = headers)
        p3_out.writeheader()
        p3_out.writerows(rows)

        return

def construct_csv_format(primers, headers): 
    rows = []
    for primer in primers:
        primers[primer]['primer'] = primer
        
        del primers[primer]['primer_start'] 
        del primers[primer]['primer_end'] 
        del primers[primer]['coords'] 

        rows.append(primers[primer])

    return rows

def construct_bed_format(primers, chromo, slice_coords, strand):
    rows = []
    for primer in primers:
        primer_data = primers[primer]

        #chr,chrStart,chrEnd,name,score,strand
        #Score unknown until iPCRess
        row = [
            chromo,
            primer_data['primer_start'],
            primer_data['primer_end'],
            primer,
            '42',
            strand
        ]
        rows.append(row)
    return rows

def export_to_bed(bed_rows):
    p3_bed = BedTool(bed_rows)
    p3_bed.saveas('./p3_output.bed')

    return
