from os import path
from collections import defaultdict
import csv
import pprint
import pdb
import re

from utils.file_system import write_to_text_file, FolderCreator
from utils.exceptions import OutputError, FolderCreatorError

class PrimerDesigner():
    def __init__(self):
        self.primer_pairs = []

    def get_primer_pairs(self):
        return self.primer_pairs

    def append_pair(self, pair):
        self.primer_pairs.append(pair)
        return

class PrimerPair():
    def __init__(self, data):
        self.pair = data['pair']
        self.score = data['score']
        self.product_size = int
        self.left = Primer(data['left'])
        self.right = Primer(data['right'])

class Primer():
    def __init__(self, primer_data):
        self.chromosome = primer_data['chromosome']
        self.chr_start = primer_data['chr_start']
        self.chr_end = primer_data['chr_end']
        self.seq = primer_data['seq']
        self.melting_temp = primer_data['melting_temp']
    
    def __getitem__(self, item):
        return getattr(self, item)


def transform_primer_pairs(primer_designer, data) -> PrimerDesigner():
    pairs = iterate_design(data)
    build_pair_classes(primer_designer, pairs)

    return primer_designer

def build_pair_classes(targeton_designer, pairs):
    for pair in pairs:
        left = extract_primer_data(pairs[pair]['F'])
        right = extract_primer_data(pairs[pair]['R'])

        pair_class = PrimerPair({
            'pair' : pair,
            'score' : 0,
            'left' : Primer(left),
            'right' : Primer(right),
        })

        targeton_designer.append_pair(pair_class)
    
    return

def extract_primer_data(data):
    record = {
        'chromosome' : data['chrom'],
        'chr_start' : data['primer_start'],
        'chr_end' : data['primer_end'],
        'seq' : data['sequence'],
        'melting_temp' : data['tm'],
    }
    
    del data['chrom']

    return record
        
def iterate_design(data):
    pairs = defaultdict(dict)

    for pair in data:
        if (pair['primers']):
            map_primer_data(pairs, pair['primers'], pair['chrom'])
    
    return pairs

def map_primer_data(pairs, data, chrom):
    for primer in data:
        #ENSE00000769557_HG8_16_LibAmpR_0
        #(ENSE00000769557_HG8_16_LibAmp)(R)_(0)
        match = re.search("^(\w+LibAmp)([F|R])\_(\d+)$", primer)
        pair_key = match.group(1) + '_' + match.group(3)
        side = match.group(2)
        pairs[pair_key][side] = data[primer]
        pairs[pair_key][side]['chrom'] = chrom

    return pairs
