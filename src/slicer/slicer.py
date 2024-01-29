#!/usr/bin/env python3
import csv
import re

from pybedtools import BedTool
from pybedtools.helpers import BEDToolsError

from utils.exceptions import SlicerError


class Slicer:
    def __init__(self):
        pass

    def get_slices(self, params):
        try:
            input_bed = params['bed']
            if params['1b']:
                input_bed = self.handle_one_based_input(params['bed'])
            bed = BedTool(input_bed)

            slice_bed = BedTool(self.get_slice_data(bed, params))
            # return named, coords slice sequences on specified strand

            return self.get_seq(slice_bed, params['fasta'])

        except Exception as err:
            raise SlicerError(f'Unexpected error occurred: {err}')

    def handle_one_based_input(self, input_bed):
        adjusted_tsv = []
        with open(input_bed) as file:
            tsv = csv.reader(file, delimiter="\t")
            adjusted_tsv = self.decrement_one_based_starts(tsv, adjusted_tsv)
        return adjusted_tsv

    def get_slice_data(self, bed, params):
        slices = []
        count = 1
        for exon in bed:
            name = exon.name if exon.name != '.' else f'region{count}'
            slices.extend(self._generate_slice_data(exon, name, params))
            count += 1
        return slices

    @staticmethod
    def get_seq(slice_bed, fasta_param):
        seq_options = {"fi": fasta_param, "s": True, "name+": True}

        try:
            seq = slice_bed.sequence(**seq_options)
        except BEDToolsError as bed_err:
            if not re.search(r'\*{5}ERROR:\ Unrecognized parameter: -name\+\ \*{5}', bed_err.args[1]):
                template = "PyBEDTools exited with err type {0}. Arguments:\n{1!r}"
                message = template.format(type(bed_err).__name__, bed_err.args[1])
                raise BEDToolsError(bed_err, message)
            del seq_options['name+']
            seq_options['name'] = True
            seq = slice_bed.sequence(**seq_options)

        return seq

    @staticmethod
    def decrement_one_based_starts(tsv, new_tsv):
        # BED is only 0-based on the start thus only need to edit column 1
        for row in tsv:
            row[1] = str(int(row[1]) - 1)
            new_tsv.append(row)
        return new_tsv

    @staticmethod
    def _generate_slice_data(exon, exon_name, params):
        slices = []
        start = exon.start - params['flank_5']
        end = start + params['length']
        count = 1
        while end <= (exon.end + params['flank_3']):
            slice_name = f'{exon_name}_{count}'
            slices.append((exon.chrom, start, end, slice_name, exon.score, exon.strand))
            start += params['offset']
            end += params['offset']
            count += 1
        return slices
