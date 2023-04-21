SRA_TOOL_KIT_URL = 'https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz'

BENCHMARKING_CONFIG = {
    'workdir': '/data/benchmarking',
    'sample_info': {
        'input': {
            'genome_reference': 'hg38',
            'patient_id': 'HCC_1395',
            'samples': [
                {'ID': 'SRR7890852', 'name': 'WES_NV_T_1', 'origin': 'tumor'},
                {'ID': 'SRR7890847', 'name': 'WES_NV_N_1', 'origin': 'normal'},
            ]
        },
        'SEQC2': {
            'reference': {
                'url': 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest',
                'files': [
                    'High-Confidence_Regions_v1.2.bed',
                    'high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz', 'high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz',
                    'high-confidence_sINDEL_in_HC_regions_v1.2.vcf.gz', 'high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz',
                    'sINDEL.MDKT.superSet.v1.2.vcf.gz',
                    'sSNV.MSDUKT.superSet.v1.2.vcf.gz',
                    'README.md'
                ]
            },
        }
    },
    'SEQC2_reference_variant_calls_dir': 'reference_variant_calls',
    'neoantigen_discovery_pipelines': {
        'old': {
            'workdir': 'old',
            'image': 'neoantigen-vaccine-pipeline:latest'
        },
        'new': {
            'workdir': 'old',
            'image': 'neoantigen-vaccine-pipeline-new:latest'
        },
    },
    'variant_call_vcf_copies': {
        'workdir':'variant_calls',
        'subdir': ['old','new']
    },
    'final_output_dir': 'final_output',
}