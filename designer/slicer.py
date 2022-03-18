import pandas as pd

def import_bed_file(file):
    col_names = ['chr', 'start', 'end']
    return pd.read_csv(file, delimiter='\t', names=col_names)
