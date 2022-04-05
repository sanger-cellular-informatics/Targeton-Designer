import primer3
import json
import pdb
import re
import collections

def primer3_runner():
    design = primer3_design()
    primers = locate_primers(design)
    return primers

def primer3_design():
    design = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG',
        'SEQUENCE_INCLUDED_REGION': [36,342]
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],
    })
    return design

def locate_primers(design):
    primer_keys = design.keys()
    primers = collections.defaultdict(dict)
    for key in primer_keys:
        match = re.search(r'^(PRIMER_(LEFT|RIGHT)_\d+)(\_(\S+))$', key)
        if match:
            primer_name = match.group(1).lower()
            primer_field = match.group(4).lower()
            primers[primer_name][primer_field] = design[key]
    return primers
