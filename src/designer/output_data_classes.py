from dataclasses import dataclass


@dataclass
class OutputFilesData:
    dir: str

    def fields(self):
        return list(self.__dataclass_fields__.keys())


@dataclass
class SlicerOutputData(OutputFilesData):
    bed: str = ''
    fasta: str = ''


@dataclass
class PrimerOutputData(OutputFilesData):
    bed: str = ''
    csv: str = ''
    discarded_csv: str = ''
    optimal_primer_pairs_csv: str = ''
    primer_pairs_for_ipcress: str = ''


@dataclass
class PrimerDesignerOutputData(OutputFilesData):
    csv: str = ''
    json: str = ''


@dataclass
class TargetonCSVData(OutputFilesData):
    csv: str = ''


@dataclass
class ScoringOutputData(OutputFilesData):
    tsv: str = ''


@dataclass
class DesignOutputData(OutputFilesData):
    p3_bed: str = ''
    p3_csv: str = ''
    pd_csv: str = ''
    pd_json: str = ''
    targeton_csv: str = ''
    scoring_tsv: str = ''
    