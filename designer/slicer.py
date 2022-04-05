import pandas as pd

def import_bed_file(file):
    col_names = ['chr', 'start', 'end']
    return pd.read_csv(file, delimiter='\t', names=col_names)

def _generate_slice_coordinates(exon, flank_5, flank_3, length, offset):
    slices = {}
    targeton_start =  exon['start'] - flank_5
    targeton_end = targeton_start + length
    while targeton_end <= (exon['end'] + flank_3):
        slices[targeton_start] = targeton_end
        targeton_start += offset
        targeton_end += offset
    return slices

def add_slice_coordinates(exons, flank_5, flank_3, length, offset):
    exons['slices'] = exons.apply(
        _generate_slice_coordinates,
        axis=1,
        args=(flank_5, flank_3, length, offset)
    )
    return exons
