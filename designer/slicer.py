from pybedtools import BedTool

def _generate_slice_coordinates(exon, params):
    slices = []
    start = exon.start - params['flank_5']
    end = start + params['length']
    while end <= (exon.end + params['flank_3']):
        slices.append((exon.chrom, start, end))
        start += params['offset']
        end += params['offset']
    return slices

def get_slice_coordinates(bed, params):
    slices = []
    for exon in bed:
        slices.extend(_generate_slice_coordinates(exon, params))
    return slices

def get_slice_sequences(bed, fasta):
    seqs = {}
    results = bed.sequence(fi=fasta, tab=True).print_sequence().strip()
    for row in results.split('\n'):
        name, sequence = row.split('\t')
        seqs[name] = sequence
    return seqs
