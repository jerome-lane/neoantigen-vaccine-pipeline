import os
from config import BENCHMARKING_CONFIG

WORKDIR = BENCHMARKING_CONFIG['workdir']

def create_reference_vcf_file(WORKDIR, BENCHMARKING_CONFIG):
    snvs = os.path.join(
        WORKDIR, BENCHMARKING_CONFIG['reference_variant_calls'],
        'high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz'
    )
    indels = os.path.join(
        WORKDIR, BENCHMARKING_CONFIG['reference_variant_calls'],
        'high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz'
    )

    reference_vcf_file = os.path.join(WORKDIR, "reference.vcf")
    os.system(
        f"bcftools concat --allow-overlaps {snvs} {indels} > {reference_vcf_file}")
    return reference_vcf_file
