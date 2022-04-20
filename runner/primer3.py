import primer3
import json
import pdb
import re
import os
import collections

def primer3_runner():
    design = primer3_design()
    primers = locate_primers(design)
    return primers

def primer3_design():
    p3_config_loc = os.environ.get('PRIMER3_CONFIG')
    print(p3_config_loc,flush=True)
    primer3_config = {}
    with open(os.path.join(os.path.dirname(__file__), p3_config_loc), "r+") as p3:
        primer3_config = json.load(p3)

    design = primer3.bindings.designPrimers({
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG',
        'SEQUENCE_INCLUDED_REGION': [36,342]
    }, primer3_config)

    return design

def locate_primers(design):
    primer_keys = design.keys()
    primers = collections.defaultdict(dict)
    for key in primer_keys:
        match = re.search(r'^(PRIMER_(LEFT|RIGHT)_\d+)\_(\S+)$', key)
        if match:
            primer_name = match.group(1).lower()
            primer_field = match.group(3).lower()
            primers[primer_name][primer_field] = design[key]
    return primers
