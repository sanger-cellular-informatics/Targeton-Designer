import pysam


def check_variant_in_region(chrom: str, start: int, end: int) -> bool:

    with pysam.VariantFile('examples/HAP1_filter_examples.vcf', 'r') as vcf_file:
        for variant in vcf_file:
            if variant.chrom == chrom and start <= variant.pos <= end:
                return True

    vcf_file.close()
    return False
