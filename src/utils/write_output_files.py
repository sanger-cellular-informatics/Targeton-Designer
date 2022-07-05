from os import path

from .file_system import FolderCreator

class OutputError(Exception):
    pass

def timestamped_dir(prefix):
    try:
        FolderCreator.create_timestamped(prefix)
    except FolderCreatorError as err:
        raise OutputError(f'Error creating folder: {err}')
    return FolderCreator.get_dir()

def write_slicer_output(dir_name, slices):
    BED_OUTPUT = 'slices.bed'
    FASTA_OUTPUT = 'fasta.bed'

    dir = timestamped_dir(dir_name)

    slices.saveas(path.join(dir, BED_OUTPUT))
    slices.save_seqs(path.join(dir, FASTA_OUTPUT))
    print('Slice files saved')