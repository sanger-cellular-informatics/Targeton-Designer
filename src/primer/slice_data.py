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
            'SEQUENCE_TEMPLATE': (
                self.surrounding_region
                if self.flanking_region and self.flanking_region > 0
                else self.bases
            ),
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': primer_pair_region
        }

    @property
    def surrounding_region(self) -> str:
        flanking = self.flanking_region

        return get_seq_from_ensembl_by_coords(
            chromosome=self.chromosome,
            start = max(1, self.start - flanking),
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
            # Chromosomes 1-22, X, Y and MT
            # Prefix for chromosome (chr or ch) is optional
            match = re.search(r'^(\w+)::(?:chr|ch|)([1-9]|1[0-9]|2[0-2]|X|Y|MT):(\d+)\-(\d+)\(([+-\.]{1})\)$', first_row.id)
            if not match:
                raise ValueError(f"The sequence ID '{first_row.id}' does not match the expected format.")

            slice_data = SliceData(
                name=match.group(1),
                start=int(match.group(3)),
                end=int(match.group(4)),
                strand=match.group(5),
                chromosome=match.group(2),
                bases=str(first_row.seq),
                flanking_region=flanking,
                exclusion_region=exclusion_region)

            if next(rows, None) is not None:
                logger.warning(f"The FASTA file '{fasta}' contains more than one pre-targeton. "
                               "Only the first pre-targeton is taken.")

        return slice_data

    @staticmethod
    def get_slice_from_region(targeton_id: str,
                              region: str,
                              strand: str,
                              flanking: int,
                              exclusion_region: int) -> 'SliceData':

        # region chr19:54100-541200
        match = re.search(r'^chr([1-9]|1[0-9]|2[0-2]|X|Y|MT):(\d+)\-(\d+)$', region)
        if not match:
            msg = f"region '{region}' does not match the expected format e.g. chr19:54100-54200"
            raise ValueError(msg)

        chromosome = match.group(1)
        target_start = int(match.group(2))
        target_end = int(match.group(3))

        if flanking and flanking > 0:
            extended_start = max(1, target_start - flanking)
            extended_end = target_end + flanking
        else:
            extended_start = target_start
            extended_end = target_end
            flanking = 0

        seq = get_seq_from_ensembl_by_coords(
            chromosome=chromosome,
            start=extended_start,
            end=extended_end,
            strand=strand
        )

        # Logging: internal vs extended region & sequences
        logger.info(
            f"Target region (internal): chr{chromosome}:{target_start}-{target_end} (strand {strand})"
        )
        logger.info(
            f"Flanking region (bp): {flanking}"
        )
        logger.info(
            f"Extended region used for template: chr{chromosome}:{extended_start}-{extended_end}"
        )
        logger.info(
            f"Internal sequence length: {target_end - target_start + 1}"
        )
        logger.info(
            f"Extended sequence length: {len(seq)}"
        )

        slice_data = SliceData(
            name=targeton_id,
            start=extended_start,
            end=extended_end,
            strand=strand,
            chromosome=chromosome,
            bases=seq,
            flanking_region=flanking,
            exclusion_region=exclusion_region
        )

        return slice_data