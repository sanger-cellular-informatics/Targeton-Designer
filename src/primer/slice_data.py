import re
from Bio import SeqIO

from primer.ensembl import get_seq_from_ensembl_by_coords

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class SliceData:
    def __init__(self, name: str, start: str, end: str, strand: str, chrom: str, bases: str):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.chrom = chrom
        self.bases = bases
        self.targeton_id = name[0:4]

    def __repr__(self):
        return (f'SliceData({self.name}, {self.targeton_id}, {self.start}, {self.end}, {self.strand}, {self.chrom},'
                f' {self.bases})')

    @property
    def p3_input(self):
        return {
            'SEQUENCE_ID': self.name,
            'SEQUENCE_TEMPLATE': self.bases,
        }

    @property
    def surrounding_region(self) -> str:
        surrounding_band = 1000

        return get_seq_from_ensembl_by_coords(
            chromosome=self.chrom,
            start=int(self.start) - surrounding_band,
            end=int(self.end) + surrounding_band
        )

    @staticmethod
    def get_first_slice_data(fasta: str) -> 'SliceData':
        with open(fasta) as fasta_data:
            rows = SeqIO.parse(fasta_data, 'fasta')

            first_row = next(rows, None)
            if first_row is None:
                raise ValueError(f"The FASTA file '{fasta}' is empty and no contains expected sequence pattern.")

            # Name::Chr:Start-End(Strand)
            # ENSE00000769557_HG8_1::1:42929543-42929753
            match = re.search(r'^(\w+)::((chr)?\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', first_row.id)

            if not match:
                raise ValueError(f"The sequence ID '{first_row.id}' does not match the expected format.")

            slice_data = SliceData(
                name=match.group(1),
                start=match.group(4),
                end=match.group(5),
                strand=match.group(6),
                chrom=match.group(2),
                bases=str(first_row.seq),
            )

        if next(rows, None) is not None:
            logger.warning(f"The FASTA file '{fasta}' contains more than one pre-targeton. "
                           f"Only the first pre-targeton is taken.")

        return slice_data
