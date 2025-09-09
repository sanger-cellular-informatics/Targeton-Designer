import re
from Bio import SeqIO

from primer.ensembl import get_seq_from_ensembl_by_coords

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class SliceData:
    def __init__(self, name: str, start: int, end: int, strand: str, chromosome: str, bases: str, region_padding: int):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.chromosome = chromosome
        self.bases = bases
        self.targeton_id = name[0:4]
        self.region_padding = region_padding

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
        print("surrounding", self.surrounding_region)
        return {
            'SEQUENCE_ID': self.name,
            'SEQUENCE_TEMPLATE': self.surrounding_region,
        }

    @property
    def surrounding_region(self) -> str:
        padding = self.region_padding

        return get_seq_from_ensembl_by_coords(
            chromosome=self.chromosome,
            start=self.start - padding,
            end=self.end + padding
        )
    
    @property
    def padded_start_coordinate(self) -> int:
        return self.start - self.region_padding
    
    @property
    def padded_end_coordinate(self) -> int:
        return self.end + self.region_padding

    @staticmethod
    def get_first_slice_data(fasta: str, padding: int) -> 'SliceData':
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
                region_padding=padding
            )

            if next(rows, None) is not None:
                logger.warning(f"The FASTA file '{fasta}' contains more than one pre-targeton. "
                               "Only the first pre-targeton is taken.")

        return slice_data
