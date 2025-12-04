import re
from Bio import SeqIO

from primer.ensembl import get_seq_from_ensembl_by_coords
from primer.flanking import build_flanked_slice
from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


# Name::Chr:Start-End(Strand)
# ENSE00000769557_HG8_1::1:42929543-42929753
# Chromosomes 1-22, X, Y and MT
# Prefix for chromosome (chr or ch) is optional
PRETARGETON_RE = re.compile(
    r'^(\w+)::(?:chr|ch|)?([1-9]|1[0-9]|2[0-2]|X|Y|MT):(\d+)-(\d+)\(([+-\.])\)$'
)

# Region mode: chr19:54100-54200
REGION_RE = re.compile(
    r'^chr([1-9]|1[0-9]|2[0-2]|X|Y|MT):(\d+)-(\d+)$'
)


def parse_pretargeton_header(header: str) -> dict:
    """
    Parse FASTA headers like:
      NAME::chr1:12345-67890(+)

    Returns:
        dict(name, chromosome, start, end, strand)

    Raises:
        ValueError: if the header does not match the expected format.
    """
    match = PRETARGETON_RE.match(header)
    if not match:
        raise ValueError(f"The sequence ID '{header}' does not match the expected format.")
    return {
        "name": match.group(1),
        "chromosome": match.group(2),
        "start": int(match.group(3)),
        "end": int(match.group(4)),
        "strand": match.group(5),
    }


def parse_region_string(region: str) -> dict:
    """
    Parse region strings like:
      chr19:54100-54200

    Returns:
        dict(chromosome, start, end)

    Raises:
        ValueError: if the region does not match the expected format.
    """
    match = REGION_RE.match(region)
    if not match:
        msg = f"region '{region}' does not match the expected format e.g. chr19:54100-54200"
        raise ValueError(msg)
    return {
        "chromosome": match.group(1),
        "start": int(match.group(2)),
        "end": int(match.group(3)),
    }


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
            # left_start = 0
            # left_len = primer_region_length

            # right_start = len(self.bases) - primer_region_length
            # right_len = primer_region_length

            # primer_pair_region = [
            #     left_start, left_len,
            #     right_start, right_len,
            # ]
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
    def get_first_slice_data(fasta: str,
                             flanking: int,
                             exclusion_region:int) -> 'SliceData':

        with open(fasta) as fasta_data:
            rows = SeqIO.parse(fasta_data, 'fasta')
            first_row = next(rows, None)
            if first_row is None:
                raise ValueError(f"Unable to parse the FASTA file '{fasta}'")

            fields = parse_pretargeton_header(first_row.id)

            name = fields["name"]
            chromosome = fields["chromosome"]
            target_region_start = fields["start"]
            target_region_end = fields["end"]
            strand = fields["strand"]
            fasta_seq = str(first_row.seq)

            if next(rows, None) is not None:
                logger.warning(f"The FASTA file '{fasta}' contains more than one pre-targeton. "
                               "Only the first pre-targeton is taken.")

        if flanking == 0:
            return SliceData(
                name=name,
                start=target_region_start,
                end=target_region_end,
                strand=strand,
                chromosome=chromosome,
                bases=fasta_seq,
                flanking_region=flanking,
                exclusion_region=exclusion_region)
        
        # flanking > 0
        flanked_start, flanked_end, flanked_seq = build_flanked_slice(
            chromosome=chromosome,
            target_region_start=target_region_start,
            target_region_end=target_region_end,
            strand=strand,
            flanking=flanking,
            mode="FASTA",
            note="original FASTA template discarded; sequence retrieved from Ensembl",
        )

        return SliceData(
            name=name,
            start=flanked_start,
            end=flanked_end,
            strand=strand,
            chromosome=chromosome,
            bases=flanked_seq,
            flanking_region=flanking,
            exclusion_region=exclusion_region
        )

    @staticmethod
    def get_slice_from_region(targeton_id: str,
                              region: str,
                              strand: str,
                              flanking: int,
                              exclusion_region: int) -> 'SliceData':

        fields = parse_region_string(region)

        chromosome = fields["chromosome"]
        target_region_start = fields["start"]
        target_region_end = fields["end"]
        
        flanked_start, flanked_end, flanked_seq = build_flanked_slice(
            chromosome=chromosome,
            target_region_start=target_region_start,
            target_region_end=target_region_end,
            strand=strand,
            flanking=flanking,
            mode="REGION",
        )

        return SliceData(
            name=targeton_id,
            start=flanked_start,
            end=flanked_end,
            strand=strand,
            chromosome=chromosome,
            bases=flanked_seq,
            flanking_region=flanking if flanking > 0 else 0,
            exclusion_region=exclusion_region
        )
