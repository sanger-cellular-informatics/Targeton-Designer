import primer3
import json
import re
import os
import collections
import csv
import sys
import argparse

from pybedtools import BedTool
from Bio import SeqIO

def parse_args(args):
    parser = argparse.ArgumentParser(
        description='Slice analysis using Primer3')
    parser.add_argument('--seq',
        help='FASTA file from the slicer tool containing seqs and IDs')
    parser.add_argument('--bed',
        help='BED file from the slicer tool containing coords, strand and IDs')
    parser.add_argument('--ref',
        help='Genomic Reference file')
    parser.add_argument('--dir',
        help='Output folder location')
    return parser.parse_args(args)

def primer3_runner(params):
    print('Reading FA file')
    design_inputs = read_input_fasta(params)
    print('Designing primers for the region')
    designs = primer3_design(design_inputs)
    print('Naming primers')
    slices = locate_primers(designs)
    print('Exporting to BED and CSV')
    export_primers(slices, params['dir'])
 
    return slices

def read_input_fasta(params):
    rows = SeqIO.parse(open(params['seq']), 'fasta')    
   
    slices = [] 
    for row in rows:
        #Name::Chr:Start-End(Strand)
        #ENSE00000769557_HG8_1::1:42929543-42929753
        match = re.search(r'^(\w+)::((chr)?\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', row.id)        
        if match:
            slice_data = construct_slice_coord_dict(match)
            p3_input = {
                'SEQUENCE_ID' : slice_data['name'],
                'SEQUENCE_TEMPLATE' : str(row.seq),
            }
            slice_data['p3_input'] = p3_input
            slices.append(slice_data)

    return slices

def construct_slice_coord_dict(match):
    coord_data = {
        'name'      : match.group(1),
        'start'     : match.group(4),
        'end'       : match.group(5),
        'strand'    : match.group(6),
        'chrom'     : match.group(2),
    }
    return coord_data

def primer3_design(primer3_inputs):
    p3_config_loc = os.environ.get('PRIMER3_CONFIG')
   
    designs = [] 
    for slice_data in primer3_inputs:
        primer3_input = slice_data['p3_input']
        primer3_config = {}
        with open(os.path.join(os.path.dirname(__file__), p3_config_loc), "r") as p3:
            primer3_config = json.load(p3)

        design = primer3.bindings.designPrimers(primer3_input, primer3_config)
        slice_data['design'] = design
        designs.append(slice_data)

    return designs

def locate_primers(designs):
    slice_designs = []
    for slice_data in designs:
        design = slice_data['design']
        primer_keys = design.keys()
        primers = collections.defaultdict(dict)

        for key in primer_keys:
            primer_details = capture_primer_details(key)
            if primer_details:
                primer_id = primer_details['id']
                primer_field = primer_details['field']
                pair_number = primer_details['pair']
                
                libamp_name = name_primers(primer_details, slice_data['strand'])
                
                primer_name = slice_data['name'] + "_" + libamp_name + "_" + pair_number
                primers[primer_name][primer_field] = design[key]
                
                if primer_field == 'coords':
                    primer_coords = calculate_primer_coords(primer_details['side'], design[key], slice_data['start'])
                    primers[primer_name]['primer_start'] = primer_coords[0] 
                    primers[primer_name]['primer_end'] = primer_coords[1] 
        del slice_data['design']
        slice_data['primers'] = primers
        slice_designs.append(slice_data)
    return slice_designs

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
        '+' : fwd_primers,
        '-' : rev_primers,
    }

    primer_name = names[strand][primer_details['side']]

    return primer_name

def export_primers(slices, output_dir):
    bed_rows = construct_bed_format(slices)
    export_to_bed(bed_rows, output_dir)
    export_to_csv(slices, output_dir)

def export_to_csv(slices, output_dir):
    headers = ['primer', 'sequence', 'tm', 'gc_percent', 'penalty', 'self_any_th', 'self_end_th', 'hairpin_th', 'end_stability']
    rows = construct_csv_format(slices, headers)

    path = output_dir + '/p3_output.csv'
    with open(path, "w") as p3_fh:
        p3_out = csv.DictWriter(p3_fh, fieldnames = headers)
        p3_out.writeheader()
        p3_out.writerows(rows)

        return

def construct_csv_format(slices, headers): 
    rows = []

    for slice_data in slices:
        primers = slice_data['primers']
        for primer in primers:
            primers[primer]['primer'] = primer
            
            del primers[primer]['primer_start'] 
            del primers[primer]['primer_end'] 
            del primers[primer]['coords'] 

            rows.append(primers[primer])

    return rows

def construct_bed_format(slices):
    rows = []
    for slice_data in slices:
        primers = slice_data['primers']
        for primer in primers: 
            primer_data = primers[primer]
            #chr,chrStart,chrEnd,name,score,strand
            #Score unknown until iPCRess
            row = [
                slice_data['chrom'],
                primer_data['primer_start'],
                primer_data['primer_end'],
                primer,
                '0',
                slice_data['strand']
            ]
            rows.append(row)
    return rows

def export_to_bed(bed_rows, output_dir):
    p3_bed = BedTool(bed_rows)
    path = output_dir + '/p3_output.bed'
    p3_bed.saveas(path)

    return

def main(params):
    if os.environ.get("PRIMER3_CONFIG") is None:
        os.environ["PRIMER3_CONFIG"] = "./primer3_config.json"
    primer3_runner(params)

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    main(vars(args))

