import pysam


def check_variant_in_region(chromosome: str, start: int, end: int) -> bool:

    with pysam.VariantFile('examples/HAP1_filter_examples.vcf', 'r') as vcf_file:
        for variant in vcf_file:
            if variant.chrom == chromosome and start <= variant.pos <= end:
                return True

    return False
