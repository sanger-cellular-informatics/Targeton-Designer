import re
from typing import List
from Bio import SeqIO


class SliceData:
    def __init__(self, name: str, start: str, end: str, strand: str, chrom: str, bases: str):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.chrom = chrom
        self.bases = bases

    @property
    def p3_input(self):
        return {
            'SEQUENCE_ID': self.name,
            'SEQUENCE_TEMPLATE': self.bases,
        }


def parse_slice(match) -> dict:
    coord_data = {
        'name': match.group(1),
        'start': match.group(5),
        'end': match.group(6),
        'strand': match.group(7),
        'chrom': match.group(3),
    }
    return coord_data

def parse_fasta(fasta: str) -> List[SliceData]:
    with open(fasta) as fasta_data:
        rows = SeqIO.parse(fasta_data, 'fasta')

        slices = []
        for row in rows:
            # Name::Chr:Start-End(Strand)
            # ENSE00000769557_HG8_1::1:42929543-42929753
            match = re.search(
              r'^(ENSE\d+_HG\d+)(_\d+)::((chr)?\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', row.id)
            if match:
                parsed_id = parse_slice(match)

                slice = SliceData(
                    parsed_id["name"],
                    parsed_id["start"],
                    parsed_id["end"],
                    parsed_id["strand"],
                    parsed_id["chrom"],
                    str(row.seq),
                )

                slices.append(slice)

    return slices