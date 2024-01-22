import re
from Bio import SeqIO


def construct_slice_coord_dict(match) -> dict:
    coord_data = {
        'name': match.group(1),
        'start': match.group(4),
        'end': match.group(5),
        'strand': match.group(6),
        'chrom': match.group(2),
    }
    return coord_data

def parse_fasta(fasta: str) -> str:
    with open(fasta) as fasta_data:
        rows = SeqIO.parse(fasta_data, 'fasta')

        slices = []
        for row in rows:
            # Name::Chr:Start-End(Strand)
            # ENSE00000769557_HG8_1::1:42929543-42929753
            match = re.search(
                r'^(\w+)::((chr)?\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', row.id)
            if match:
                slice_data = construct_slice_coord_dict(match)
                p3_input = {
                    'SEQUENCE_ID': slice_data['name'],
                    'SEQUENCE_TEMPLATE': str(row.seq),
                }
                slice_data['p3_input'] = p3_input
                slices.append(slice_data)

    return slices