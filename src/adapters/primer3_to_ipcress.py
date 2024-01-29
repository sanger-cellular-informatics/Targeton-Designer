import os
import re
from collections import defaultdict

from utils.file_system import read_csv_to_list_dict


class Primer3ToIpcressAdapter:
    def __init__(self):
        self.primer_data = []
        self.formatted_primers = []

    def prepare_input(self, csv, min, max, dir):
        primer_data = self.extract_primer_data(csv)
        self.primer_data = primer_data

        formatted_primers = self.get_ipcress_input(primer_data, min, max, dir)
        self.formatted_primers = formatted_primers

        return formatted_primers

    def get_ipcress_input(self, primer_data, min, max, dir):
        formatted_primers = self.format_ipcress_primers(min, max, primer_data)

        return formatted_primers

    @staticmethod
    def extract_primer_data(p3_csv):
        csv_obj = read_csv_to_list_dict(p3_csv)
        primer_data = defaultdict(lambda: defaultdict(dict))

        for row in csv_obj:
            # Capture primer name and orientation
            # ENSE00003571441_HG6_6_LibAmpR_0
            match = re.search(r'^(\w+_LibAmp)(F|R)_(\d+)$', row['primer'])
            if match:
                key = match.group(1) + '_' + match.group(3)
                primer_data[key][match.group(2)]['seq'] = row['sequence']
                primer_data[key][match.group(2)]['start'] = row['primer_start']

        return primer_data

    @staticmethod
    def format_ipcress_primers(min_amp, max_amp, primers):
        ipcress_input = []
        rows = primers.keys()

        for key in rows:
            left = primers[key]['F']['seq']
            right = primers[key]['R']['seq']

            line = ' '.join([key, left, right, min_amp, max_amp])
            ipcress_input.append(line)

        return ipcress_input
