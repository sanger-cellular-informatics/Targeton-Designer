import pysam
from pathlib import Path
from Bio.Seq import Seq


def get_seq_from_fasta_by_coords(
    fasta_file: str | Path,
    chromosome: str,
    start: int,
    end: int,
    strand: str = "+",
) -> str:
    """Fetch a 1-based, inclusive sequence range from an indexed FASTA file."""

    if start < 1 or end < start:
        raise ValueError("Coordinates must be 1-based, with end >= start.")

    if strand not in {"+", "-"}:
        raise ValueError("strand must be '+' or '-'.")

    fasta_path = Path(fasta_file)

    try:
        with pysam.FastaFile(str(fasta_path)) as fasta:
            # pysam uses 0-based, end-exclusive coordinates.
            sequence = fasta.fetch(chromosome, start - 1, end).upper()
    except (OSError, ValueError, KeyError) as exc:
        raise ValueError(
            f"Could not fetch {chromosome}:{start}-{end} from {fasta_path}. "
            "Check the FASTA index, chromosome name, and coordinates."
        ) from exc

    if strand == "-":
        # reverse complement
        sequence = str(Seq(sequence).reverse_complement())

    return sequence
