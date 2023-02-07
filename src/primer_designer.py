from os import path
from collections import defaultdict
import re

class PrimerDesigner():
    def __init__(self):
        self.primer_pairs = []

    def get_primer_pairs(self):
        return self.primer_pairs

    def append_pair(self, pair):
        self.primer_pairs.append(pair)
        
    def get_fields(self):
        fields = []
        for pair in self.primer_pairs:
            fields.extend(list(pair.get_paired_dict().keys()))
        fields = list(set(fields))
        return fields

class PrimerPair():
    def __init__(self, data):
        self.pair = data['pair']
        self.score = data['score']
        self.product_size = int
        self.left = Primer(data['left'])
        self.right = Primer(data['right'])

    def get_paired_dict(self):
        return vars(self)
    
    def get_fields(self):
        return list(self.attrs().keys())

class Primer():
    def __init__(self, primer_data):
        self.chromosome = primer_data['chromosome']
        self.chr_start = primer_data['chr_start']
        self.chr_end = primer_data['chr_end']
        self.seq = primer_data['seq']
        self.melting_temp = primer_data['melting_temp']
    
    def __getitem__(self, item):
        return getattr(self, item)
    
    def get_fields(self):
        return list(self.attrs().keys())


def prepare_primer_designer(primer_designer: PrimerDesigner, primers: list) -> PrimerDesigner:
    pairs = iterate_design(primers)
    primer_designer = build_pair_classes(primer_designer, pairs)

    return primer_designer

def build_pair_classes(primer_designer:PrimerDesigner, pairs: defaultdict) -> PrimerDesigner:
    for pair in pairs:
        left = extract_primer_data(pairs[pair]['F'])
        right = extract_primer_data(pairs[pair]['R'])

        pair_class = PrimerPair({
            'pair' : pair,
            'score' : 0,
            'left' : Primer(left),
            'right' : Primer(right),
        })

        primer_designer.append_pair(pair_class)
    
    return primer_designer

def extract_primer_data(data: dict) -> dict:
    record = {
        'chromosome' : data['chrom'],
        'chr_start' : data['primer_start'],
        'chr_end' : data['primer_end'],
        'seq' : data['sequence'],
        'melting_temp' : data['tm'],
    }
    
    # For ipcress compatibility.
    del data['chrom']

    return record
        
def iterate_design(primers: list) -> defaultdict:
    pairs = defaultdict(dict)

    for pair in primers:
        if (pair['primers']):
            pairs.update(map_primer_data(pairs, pair['primers'], pair['chrom']))
    
    return pairs

def map_primer_data(pairs: defaultdict, primers: defaultdict, chrom: str):
    for primer in primers:
        #ENSE00000769557_HG8_16_LibAmpR_0
        #(ENSE00000769557_HG8_16_LibAmp)(R)_(0)
        match = re.search("^(\w+LibAmp)([F|R])\_(\d+)$", primer)
        pair_key = match.group(1) + '_' + match.group(3)
        side = match.group(2)
        pairs[pair_key][side] = primers[primer]
        pairs[pair_key][side]['chrom'] = chrom

    return pairs
