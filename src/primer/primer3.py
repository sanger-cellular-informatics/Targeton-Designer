import primer3
import json
import re
import os
import collections

from Bio import SeqIO
from Bio.Seq import Seq

from src.utils.exceptions import Primer3Error


def primer3_runner(fasta=''):
    print('Reading FA file')
    design_inputs = read_input_fasta(fasta)
    print('Designing primers for the region')
    designs = primer3_design(design_inputs)
    print('Naming primers')
    slices = locate_primers(designs)

    return slices


def read_input_fasta(fasta):
    rows = SeqIO.parse(open(fasta), 'fasta')

    slices = []
    for row in rows:
        # Name::Chr:Start-End(Strand)
        # ENSE00000769557_HG8_1::1:42929543-42929753
        match = re.search(r'^(\w+)::((chr)?\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', row.id)
        if match:
            slice_data = construct_slice_coord_dict(match)
            p3_input = {
                'SEQUENCE_ID': slice_data['name'],
                'SEQUENCE_TEMPLATE': str(row.seq),
            }
            slice_data['p3_input'] = p3_input
            slices.append(slice_data)

    return slices


def construct_slice_coord_dict(match):
    coord_data = {
        'name'  : match.group(1),
        'start' : match.group(4),
        'end'   : match.group(5),
        'strand': match.group(6),
        'chrom' : match.group(2),
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
                primer = primers[primer_name]

                primer[primer_field] = design[key]
                primer['side'] = primer_details['side']

                if primer_field == 'coords':
                    primer_coords = calculate_primer_coords(primer_details['side'], design[key], slice_data['start'])
                    primer['primer_start'] = primer_coords[0]
                    primer['primer_end'] = primer_coords[1]
                    primer['strand'] = determine_primer_strands(primer_details['side'], slice_data['strand'])
                    primer['sequence'] = revcom_reverse_primer(primer['sequence'], primer['strand'])
                primers[primer_name] = primer
        del slice_data['design']
        slice_data['primers'] = primers
        slice_designs.append(slice_data)
    return slice_designs


def revcom_reverse_primer(seq, strand):
    seq_obj = Seq(seq)

    if strand == '-':
        seq_obj = seq_obj.reverse_complement()

    return seq_obj


def determine_primer_strands(side, slice_strand):
    positive = {
        'left': '+',
        'right': '-',
    }

    negative = {
        'left': '-',
        'right': '+',
    }

    strands = {
        '+': positive,
        '-': negative,
    }

    return strands[slice_strand][side]


def calculate_primer_coords(side, coords, slice_start):
    slice_start = int(slice_start)
    left_flank = {
        'start': slice_start,
        'end': slice_start + int(coords[1])
    }

    slice_end = slice_start + int(coords[0])
    right_flank = {
        'start': slice_end - int(coords[1]),
        'end': slice_end,
    }

    slice_coords = {
        'left': left_flank,
        'right': right_flank
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
            'id'    : primer_id,
            'side'  : primer_side,
            'field' : primer_field,
            'pair'  : pair_number
        }

    return result


def name_primers(primer_details, strand):
    fwd_primers = {
        'left': 'LibAmpF',
        'right': 'LibAmpR',
    }
    rev_primers = {
        'left': 'LibAmpR',
        'right': 'LibAmpF',
    }
    names = {
        '+': fwd_primers,
        '-': rev_primers,
    }

    primer_name = names[strand][primer_details['side']]

    return primer_name


def main(fasta):
    if os.environ.get("PRIMER3_CONFIG") is None:
        os.environ["PRIMER3_CONFIG"] = "./primer3_config.json"
    result = primer3_runner(fasta=fasta)

    return result


if __name__ == '__main__':
    main()
