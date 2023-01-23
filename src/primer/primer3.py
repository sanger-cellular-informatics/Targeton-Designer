import primer3
import re
import os
import collections

from Bio import SeqIO
from Bio.Seq import Seq

from utils.exceptions import Primer3Error, InvalidConfigError
from utils.file_system import parse_json


class Primer:
    def __init__(self):
        pass

        self.user_config = './config/primer3.config.json'
        self.default_config = './src/primer/primer3.config.json'

    def get_primers(self, fasta):
        config = self.get_config_data(self.default_config, self.user_config)

        result = self.primer3_runner(fasta=fasta, config=config)

        return result

    def primer3_runner(self, fasta: str, config: dict):
        print('Reading FA file')
        design_inputs = self.read_input_fasta(fasta)
        print('Designing primers for the region')
        designs = self.primer3_design(design_inputs, config)
        print('Naming primers')
        slices = self.locate_primers(designs)

        return slices

    def read_input_fasta(self, fasta):
        rows = SeqIO.parse(open(fasta), 'fasta')

        slices = []
        for row in rows:
            # Name::Chr:Start-End(Strand)
            # ENSE00000769557_HG8_1::1:42929543-42929753
            match = re.search(r'^(\w+)::((chr)?\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', row.id)
            if match:
                slice_data = self.construct_slice_coord_dict(match)
                p3_input = {
                    'SEQUENCE_ID': slice_data['name'],
                    'SEQUENCE_TEMPLATE': str(row.seq),
                }
                slice_data['p3_input'] = p3_input
                slices.append(slice_data)

        return slices

    def primer3_design(self, primer3_inputs, primer3_config):
        designs = []
        for slice_data in primer3_inputs:
            primer3_input = slice_data['p3_input']

            design = primer3.bindings.designPrimers(primer3_input, primer3_config)
            slice_data['design'] = design
            designs.append(slice_data)

        return designs

    def locate_primers(self, designs):
        slice_designs = []

        for slice_data in designs:
            design = slice_data['design']
            primer_keys = design.keys()

            primers = self.build_primers_dict(design, primer_keys, slice_data)

            del slice_data['design']
            slice_data['primers'] = primers
            slice_designs.append(slice_data)

        return slice_designs

    def build_primers_dict(self, design, primer_keys, slice_data):
        primers = collections.defaultdict(dict)

        for key in primer_keys:
            primer_details = self.capture_primer_details(key)

            if primer_details:
                libamp_name = self.name_primers(primer_details, slice_data['strand'])
                primer_name = slice_data['name'] + "_" + libamp_name + "_" + primer_details['pair']

                primers[primer_name] = self.build_primer_loci(
                    primers[primer_name], key, design,
                    primer_details, slice_data)

        return primers

    def build_primer_loci(
        self, primer, key, design, primer_details, slice_data
    ):
        
        primer_field = primer_details['field']

        primer[primer_field] = design[key]
        primer['side'] = primer_details['side']

        if primer_field == 'coords':
            primer_coords = self.calculate_primer_coords(
                primer_details['side'], design[key], slice_data['start'])

            primer['primer_start'] = primer_coords[0]
            primer['primer_end'] = primer_coords[1]
            primer['strand'] = self.determine_primer_strands(
                primer_details['side'], slice_data['strand']
            )

        return primer

    def get_config_data(self, default_config: str, user_config: str) -> dict:
        config_file = self.get_config_file(default_config, user_config)

        try:
            config_data = parse_json(config_file)

        except Exception as err:
            raise InvalidConfigError(f'Primer3 config file is not a correct JSON')

        return config_data

    @staticmethod
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

    @staticmethod
    def capture_primer_details(primer_name):
        match = re.search(r'^(primer_(left|right)_(\d+))(\_(\S+))?$', primer_name.lower())
        result = {}
        if match:
            primer_id = match.group(1)
            primer_side = match.group(2)
            pair_number = match.group(3)
            primer_field = match.group(5)
            if primer_field is None:
                primer_field = 'coords'
            result = {
                'id'    : primer_id,
                'side'  : primer_side,
                'field' : primer_field,
                'pair'  : pair_number
            }

        return result

    @staticmethod
    def construct_slice_coord_dict(match):
        coord_data = {
            'name'  : match.group(1),
            'start' : match.group(4),
            'end'   : match.group(5),
            'strand': match.group(6),
            'chrom' : match.group(2),
        }
        return coord_data

    @staticmethod
    def calculate_primer_coords(side, coords, slice_start):
        slice_start = int(slice_start)
        left_flank = {
            'start': slice_start,
            'end': slice_start + int(coords[1])
        }

        slice_end = slice_start + int(coords[0])
        right_flank = {
            'start': 1 + slice_end - int(coords[1]),
            'end': 1 + slice_end,
        }

        slice_coords = {
            'left': left_flank,
            'right': right_flank
        }

        start = slice_coords[side]['start']
        end = slice_coords[side]['end']

        return start, end

    @staticmethod
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

    @staticmethod
    def get_config_file(default_config: str, user_config:str) -> str:
        config_file_path = user_config if os.path.exists(user_config) else default_config

        return config_file_path

    @staticmethod
    def revcom_reverse_primer(seq, strand):
        seq_obj = Seq(seq)

        if strand == '-':
            seq_obj = seq_obj.reverse_complement()

        return seq_obj
