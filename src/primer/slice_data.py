import json
from pathlib import Path
from typing import Optional
import re
from Bio import SeqIO

from primer.ensembl import get_seq_from_ensembl_by_coords
from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)

CHR_LENGTHS_PATH = Path(__file__).resolve().parents[1] / "chr_lengths_grch38.json"

with CHR_LENGTHS_PATH.open() as f:
    CHR_LENGTHS_GRCh38 = json.load(f)

# ---------------------------------------------------------------------------
# Flanking helpers
# ---------------------------------------------------------------------------

def _clamp_flanked_region(flanked_start_unclamped: int,
                          flanked_end_unclamped: int,
                          chr_len: Optional[int]) -> tuple[int, int,bool,bool]:
    """
    Clamp flanked region [flanked_start_unclamped, flanked_end_unclamped] to valid chromosome bounds [1, chr_len].

    Returns (flanked_start_clamped, flanked_end_clamped), logging any adjustments.
    """
    left_clamped = False
    right_clamped = False

    if flanked_start_unclamped < 1:
        flanked_start_clamped = 1
        left_clamped = True
    else:
        flanked_start_clamped = flanked_start_unclamped

    flanked_end_clamped = flanked_end_unclamped
    if chr_len is not None and flanked_end_unclamped > chr_len:
        flanked_end_clamped = chr_len
        right_clamped = True

    return flanked_start_clamped, flanked_end_clamped, left_clamped, right_clamped

def _get_flanked_coordinates(
    chromosome: str,
    target_region_start: int,
    target_region_end: int,
    flanking: int,
) -> tuple[int, int]:
    """
    Compute flanked coordinates around a target region and clamp them
    to chromosome bounds if needed.

    Returns:
        flanked_start_clamped, flanked_end_clamped
    """
    assert flanking > 0

    chr_len = get_chromosome_length(chromosome)

    flanked_start_unclamped = target_region_start - flanking
    flanked_end_unclamped = target_region_end + flanking

    (
        flanked_start_clamped,
        flanked_end_clamped,
        left_clamped,
        right_clamped,
    ) = _clamp_flanked_region(
        flanked_start_unclamped=flanked_start_unclamped,
        flanked_end_unclamped=flanked_end_unclamped,
        chr_len=chr_len,
    )

    if left_clamped:
        logger.warning(
            "Flanking region extends beyond start of chromosome. "
            f"Requested flank start={flanked_start_unclamped}, clamped to 1. "
            "Left flank length will be reduced. "
        )

    if right_clamped and chr_len is not None:
        logger.warning(
            f"Flanking region expands beyond chromosome {chromosome} end "
            f"({flanked_end_unclamped} > {chr_len}). "
            "Flanked end has been clamped to the chromosome boundary. "
        )

    return flanked_start_clamped, flanked_end_clamped

def _build_flanked_slice(
    chromosome: str,
    target_region_start: int,
    target_region_end: int,
    strand: str,
    flanking: int,
) -> tuple[int, int, str]:
    """
    Compute flanked coordinates and fetch the corresponding sequence.

    Returns:
        flanked_start, flanked_end, sequence
    """
    flanked_start, flanked_end = _get_flanked_coordinates(
        chromosome=chromosome,
        target_region_start=target_region_start,
        target_region_end=target_region_end,
        flanking=flanking,
    )

    flanked_seq = get_seq_from_ensembl_by_coords(
        chromosome=chromosome,
        start=flanked_start,
        end=flanked_end,
        strand=strand,
    )

    return flanked_start, flanked_end, flanked_seq

def _log_flanking_summary(
    mode: str,
    chromosome: str,
    target_region_start: int,
    target_region_end: int,
    strand: str,
    flanking: int,
    flanked_start: int,
    flanked_end: int,
    flanked_seq: str,
    note: str = "",
) -> None:

    header = {
        "FASTA": "FASTA mode with auto-flanking",
        "REGION": "Region mode",
    }.get(mode, mode)

    note_line = f"\n    note:                {note}" if note else ""

    logger.info(
        f"""{header}:
    target_region:       {chromosome}:{target_region_start}-{target_region_end} ({strand})
    internal_length:     {target_region_end - target_region_start + 1} bp
    flanking_region:     {flanking} bp
    flanked_region:      {chromosome}:{flanked_start}-{flanked_end} ({strand})
    flanked_length:      {len(flanked_seq)} bp{note_line}
    sequence_len:        {len(flanked_seq)}
    sequence:            {flanked_seq}"""
    )

def get_chromosome_length(chromosome: str) -> Optional[int]:
    length = CHR_LENGTHS_GRCh38.get(chromosome)
    if length is None:
        logger.warning(
            f"Chromosome '{chromosome}' not found in CHR_LENGTHS_GRCh38. "
            "Flanking positions will not be restricted by chromosome bounds."
        )
    return length

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
    def get_first_slice_data(fasta: str,
                             flanking: int,
                             exclusion_region:int) -> 'SliceData':

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

            name = match.group(1)
            chromosome = match.group(2)
            target_region_start = int(match.group(3))
            target_region_end = int(match.group(4))
            strand = match.group(5)
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
        flanked_start, flanked_end, flanked_seq = _build_flanked_slice(
            chromosome=chromosome,
            target_region_start=target_region_start,
            target_region_end=target_region_end,
            strand=strand,
            flanking=flanking
        )

        _log_flanking_summary(
            mode="FASTA",
            chromosome=chromosome,
            target_region_start=target_region_start,
            target_region_end=target_region_end,
            strand=strand,
            flanking=flanking,
            flanked_start=flanked_start,
            flanked_end=flanked_end,
            flanked_seq=flanked_seq,
            note="original FASTA template discarded; sequence retrieved from Ensembl"
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

        # region chr19:54100-541200
        match = re.search(r'^chr([1-9]|1[0-9]|2[0-2]|X|Y|MT):(\d+)\-(\d+)$', region)
        if not match:
            msg = f"region '{region}' does not match the expected format e.g. chr19:54100-54200"
            raise ValueError(msg)

        chromosome = match.group(1)
        target_region_start = int(match.group(2))
        target_region_end = int(match.group(3))

        flanked_start, flanked_end, flanked_seq = _build_flanked_slice(
            chromosome=chromosome,
            target_region_start=target_region_start,
            target_region_end=target_region_end,
            strand=strand,
            flanking=flanking
        )

        _log_flanking_summary(
            mode="REGION",
            chromosome=chromosome,
            target_region_start=target_region_start,
            target_region_end=target_region_end,
            strand=strand,
            flanking=flanking,
            flanked_start=flanked_start,
            flanked_end=flanked_end,
            flanked_seq=flanked_seq
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
