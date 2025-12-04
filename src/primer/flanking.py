from pathlib import Path
from typing import Optional
import json
from custom_logger.custom_logger import CustomLogger

from primer.ensembl import get_seq_from_ensembl_by_coords

# Initialize logger
logger = CustomLogger(__name__)

CHR_LENGTHS_PATH = Path(__file__).resolve().parents[1] / "utils" / "chr_lengths_grch38.json"

with CHR_LENGTHS_PATH.open() as f:
    CHR_LENGTHS_GRCh38 = json.load(f)


def build_flanked_slice(
    chromosome: str,
    target_region_start: int,
    target_region_end: int,
    strand: str,
    flanking: int,
    mode: str,
    note: str = "",
) -> tuple[int, int, str]:
    """
      1. Computes flanked coordinates
      2. Fetches the sequence for those coordinates
      3. Applies strand orientation if needed
      4. Logs a flanking summary (if mode is provided)

    Returns:
        flanked_start, flanked_end, flanked_seq
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

    if mode is not None:
        _log_flanking_summary(
            mode=mode,
            chromosome=chromosome,
            target_region_start=target_region_start,
            target_region_end=target_region_end,
            strand=strand,
            flanking=flanking,
            flanked_start=flanked_start,
            flanked_end=flanked_end,
            flanked_seq=flanked_seq,
            note=note,
        )

    return flanked_start, flanked_end, flanked_seq


def _get_flanked_coordinates(
    chromosome: str,
    target_region_start: int,
    target_region_end: int,
    flanking: int,
) -> tuple[int, int]:
    """
    Compute flanked coordinates around a target region and clamp them
    to chromosome bounds if needed.
    """
    if flanking == 0:
        return target_region_start, target_region_end

    chr_len = _get_chromosome_length(chromosome)

    flanked_start_unclamped = target_region_start - flanking
    flanked_end_unclamped = target_region_end + flanking

    (
        flanked_start_clamped,
        flanked_end_clamped,
    ) = _clamp_flanked_region(
        chromosome=chromosome,
        flanked_start_unclamped=flanked_start_unclamped,
        flanked_end_unclamped=flanked_end_unclamped,
        chr_len=chr_len,
    )

    return flanked_start_clamped, flanked_end_clamped


def _clamp_flanked_region(
    chromosome: str,
    flanked_start_unclamped: int,
    flanked_end_unclamped: int,
    chr_len: int | None,
) -> tuple[int, int, bool, bool]:
    """
    Clamp the flanked region to [1, chr_len] and return:
        flanked_start_clamped, flanked_end_clamped
    """
    left_clamped = flanked_start_unclamped < 1
    right_clamped = chr_len is not None and flanked_end_unclamped > chr_len

    if left_clamped:
        logger.warning(
            "Flanking region extends beyond start of chromosome; start clamped to 1 and left flank length shortened."
        )

    if right_clamped and chr_len is not None:
        logger.warning(
            f"Flanking region expands beyond chromosome {chromosome} end;"
            f"end clamped to chromosome boundary {chr_len}. "
        )

    flanked_start_clamped = max(flanked_start_unclamped, 1)
    flanked_end_clamped = (
        min(flanked_end_unclamped, chr_len) if chr_len is not None else flanked_end_unclamped
    )

    return flanked_start_clamped, flanked_end_clamped


def _get_chromosome_length(chromosome: str) -> Optional[int]:
    length = CHR_LENGTHS_GRCh38.get(chromosome)
    if length is None:
        logger.warning(
            f"Chromosome '{chromosome}' not found in CHR_LENGTHS_GRCh38."
            "Flanking positions will not be restricted by chromosome bounds."
        )
    return length


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
        "FASTA": "FASTA mode",
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
