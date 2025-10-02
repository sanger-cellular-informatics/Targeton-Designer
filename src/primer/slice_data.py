import re
from Bio import SeqIO
import sys

from primer.ensembl import get_seq_from_ensembl_by_coords

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class SliceData:
    def __init__(self,
                 name: str,
                 start: int,
                 end: int,
                 strand: str,
                 chromosome: str,
                 bases: str,
                 flanking_region: int,
                 exclusion_region: int):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.chromosome = chromosome
        self.bases = bases
        self.targeton_id = name[0:4]
        self.flanking_region = flanking_region
        self.exclusion_region = exclusion_region

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SliceData):
            return False
        return  self.name == other.name and \
                self.start == other.start and \
                self.end == other.end and \
                self.strand == other.strand and \
                self.chromosome == other.chromosome and \
                self.bases == other.bases

    def __hash__(self):
        return hash((self.name, self.start, self.end, self.strand, self.chromosome))

    def __repr__(self):
        return (f'SliceData({self.name}, {self.targeton_id}, {self.start}, {self.end}, {self.strand}, {self.chromosome},'
                f' {self.bases})')

    @property
    def p3_input(self):
        primer_pair_region = []
        if self.flanking_region:
            primer_region_length = self.flanking_region - self.exclusion_region
            primer_pair_region = [0, primer_region_length,
                                      len(self.bases) - primer_region_length + 1, primer_region_length - 1]

        return {
            'SEQUENCE_ID': self.name,
            'SEQUENCE_TEMPLATE': self.bases,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': primer_pair_region
        }

    # Not currently in use
    @property
    def surrounding_region(self) -> str:
        flanking = self.flanking_region

        return get_seq_from_ensembl_by_coords(
            chromosome=self.chromosome,
            start=self.start - flanking,
            end=self.end + flanking
        )

    @staticmethod
    def get_first_slice_data(fasta: str, flanking: int, exclusion_region: int) -> 'SliceData':
        with open(fasta) as fasta_data:
            rows = SeqIO.parse(fasta_data, 'fasta')
            first_row = next(rows, None)
            if first_row is None:
                raise ValueError(f"Unable to parse the FASTA file '{fasta}'")

            # Name::Chr:Start-End(Strand)
            # ENSE00000769557_HG8_1::1:42929543-42929753
            # Only matches with numerical chromosome, not X, Y, and MT.
            match = re.search(r'^(\w+)::(chr\d+|ch\d+|\d+):(\d+)\-(\d+)\(([+-\.]{1})\)$', first_row.id)

            if not match:
                raise ValueError(f"The sequence ID '{first_row.id}' does not match the expected format.")

            chromosome = ''.join(filter(str.isdigit, match.group(2)))

            slice_data = SliceData(
                name=match.group(1),
                start=int(match.group(3)),
                end=int(match.group(4)),
                strand=match.group(5),
                chromosome=chromosome,
                bases=str(first_row.seq),
                flanking_region=flanking,
                exclusion_region=exclusion_region)

            if next(rows, None) is not None:
                logger.warning(f"The FASTA file '{fasta}' contains more than one pre-targeton. "
                               "Only the first pre-targeton is taken.")

        return slice_data
