"""Images of each pipeline to compare should have been built and pushed to AWS.
"""
from sklearn.metrics import classification_report
import wget
import os
import pandas as pd
import allel
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

benchmark_config = {
    'workdir': '/data/benchmarking',
    'sample_info': {
        'input': {
            'genome_reference': 'hg38',
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
            'workdir': 'new',
            'image': 'neoantigen-vaccine-pipeline-new:latest'
        },
    },
    'final_output_dir': 'final_output',
}

sample_info = benchmark_config['sample_info']
url = sample_info['SEQC2']['reference']['url']
files = sample_info['SEQC2']['reference']['files']

WORKDIR = benchmark_config['workdir']
output_directory = os.path.join(
    WORKDIR, benchmark_config['SEQC2_reference_variant_calls_dir'])
os.mkdir(output_directory, exist_ok=True)
for file in files:
    file_url = url + file
    filename = wget.download(file_url, out=output_directory)

sra_tool_kit_url = 'https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz'
wget.download(sra_tool_kit_url, out=WORKDIR)

os.system(
    f'tar -vxzf sratoolkit.tar.gz && export PATH=$PATH:{WORKDIR}/sratoolkit.3.0.2-ubuntu64/bin && cd $PWD/sratoolkit.3.0.2-ubuntu64/bin && vdb-config -i'
)

os.chdir(WORKDIR)

for sample in sample_info['input']['samples']:
    sample_id = sample['ID']
    os.system(f'prefetch {sample_id} --max-size 420000000000')
    os.system(f'fasterq-dump {sample_id}')

os.system('sh run_neoantigen_pipeline.sh')


def extract_variant_calls_compare_reference(
        tested_vcf, reference_vcf):
    """Extract variant calls coordinates from vcf files and compare results.
    Returns the results of the comparison. 0 if the variant was not called, 1 if"""
    tested_vcf_df = allel.vcf_to_dataframe(tested_vcf)
    reference_vcf_df = allel.vcf_to_dataframe(reference_vcf)
    tested_vcf_df['CHROM'] = tested_vcf_df['CHROM'].astype(str)
    reference_vcf_df['CHROM'] = reference_vcf_df['CHROM'].astype(str)
    tested_vcf_df['POS'] = tested_vcf_df['POS'].astype(int)
    reference_vcf_df['POS'] = reference_vcf_df['POS'].astype(int)
    merged_vcf_df = tested_vcf_df.merge(
        reference_vcf_df, on=['CHROM', 'POS'], how='outer',
        suffix=('_tested', '_reference')
    )
    variant_calls_pred = []
    variant_calls_true = []
    for variant in merged_vcf_df.iterrows():
        if variant[1]['REF_tested'] != variant[1]['REF_reference']:
            variant_calls_pred.append()
            variant_calls_true.append()
        if variant[1]['ALT_tested'] != variant[1]['ALT_reference']:
            variant_calls_pred.append()
            variant_calls_true.append()
    return variant_calls_pred, variant_calls_true


def build_performance_comparision_dataframe(
        benchmark_config, WORKDIR, VARIANT_CALLERS=['strelka', 'mutect']):
    """Build a dataframe with the performance of the different pipelines"""
    performance_df = pd.DataFrame()
    reference_vcf = os.path.join(
        WORKDIR, benchmark_config['reference_variant_calls'],
        'high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz')

    for pipeline, description in benchmark_config['neoantigen_discovery_pipelines'].items():
        for variant_caller in VARIANT_CALLERS:
            tested_vcf = os.path.join(
                WORKDIR, description['workdir'], 'output',
                f'{variant_caller}.vcf'
            )

            y_pred, y_true = extract_variant_calls_compare_reference(
                tested_vcf, reference_vcf)
            target_names = ['not called', 'called']
            results_dict = classification_report(
                y_true, y_pred, target_names=target_names, dict=True)
            results_dict['pipeline'] = pipeline
            results_dict['variant_caller'] = variant_caller

            performance_df = performance_df.append(
                results_dict,
                ignore_index=True
            )
        return performance_df


def plot_windowed_variant_density(pos, window_size, filename=None, title=None):
    # setup windows
    bins = np.arange(0, pos.max(), window_size)
    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2
    # compute variant density in each window
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size

    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)
    if filename:
        plt.savefig(filename)


def create_density_plot_for_all_chromosomes(
        vcf, output_dir, filename):
    for key in vcf.keys():
        if key.endswith('/variants/POS'):
            pos = allel.SortedIndex(vcf[key])
            plot_windowed_variant_density(
                pos, window_size=100000,
                filename=os.path.join(
                    output_dir, filename
                ),
                title='Variant density'
            )


def build_variant_density_plots(WORKDIR, benchmark_config, VARIANT_CALLERS):
    output_dir = benchmark_config['final_output_dir']
    reference_vcf = os.path.join(
        WORKDIR, benchmark_config['reference_variant_calls'],
        'high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz')
    reference_vcf = allel.vcf_to_dataframe(reference_vcf)
    create_density_plot_for_all_chromosomes(
        reference_vcf, output_dir, 'reference.jpg')
    for pipeline, description in benchmark_config['neoantigen_discovery_pipelines'].items():
        for variant_caller in VARIANT_CALLERS:
            tested_vcf = os.path.join(
                WORKDIR, description['workdir'], 'output',
                f'{variant_caller}.vcf'
            )
            tested_vcf = allel.vcf_to_dataframe(tested_vcf)
            create_density_plot_for_all_chromosomes(
                tested_vcf, output_dir, f'{pipeline}_{variant_caller}.jpg')


performance_df = build_performance_comparision_dataframe(
    benchmark_config, WORKDIR)
performance_df.to_excel(
    os.path.join(
        WORKDIR, benchmark_config['final_output_dir'], 'performance.xlsx')
)
