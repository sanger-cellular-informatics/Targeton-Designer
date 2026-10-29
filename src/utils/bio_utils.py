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
    """Fetch a 1-based, inclusive sequence range from an indexed FASTA file.

    The function supports FASTA contigs named either "1" or "chr1".
    """

    if start < 1 or end < start:
        raise ValueError("Coordinates must be 1-based, with end >= start.")

    if strand not in {"+", "-"}:
        raise ValueError("strand must be '+' or '-'.")

    try:
        fasta = pysam.FastaFile(str(fasta_file))
    except OSError as exc:
        raise ValueError(f"Could not open FASTA file: {fasta_file}.") from exc

    with fasta:
        references = set(fasta.references)

        if chromosome in references:
            fasta_chromosome = chromosome
        elif f"chr{chromosome}" in references:
            fasta_chromosome = f"chr{chromosome}"
        elif chromosome.startswith("chr") and chromosome[3:] in references:
            fasta_chromosome = chromosome[3:]
        else:
            raise ValueError(
                f"Chromosome {chromosome!r} was not found in the FASTA index."
            )

        chromosome_length = fasta.get_reference_length(fasta_chromosome)

        if end > chromosome_length:
            raise ValueError(
                f"Coordinates {start}-{end} exceed the length of "
                f"{fasta_chromosome} ({chromosome_length} bases)."
            )

        try:
            # pysam uses 0-based, end-exclusive coordinates.
            sequence = fasta.fetch(fasta_chromosome, start - 1, end).upper()
        except (OSError, ValueError, KeyError) as exc:
            raise ValueError(
                f"Could not fetch {chromosome}:{start}-{end} from {fasta_file}."
            ) from exc

    if strand == "-":
        # reverse complement
        sequence = str(Seq(sequence).reverse_complement())

    return sequence
