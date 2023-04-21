import os
import allel
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from config import BENCHMARKING_CONFIG


def plot_windowed_variant_density(pos, window_size, filename=None, title=None):
    """Create a plot of variant density given base pair window size."""
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
    """Create density plot for all chromosomes present in the vcf file."""
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


def build_variant_density_plots(WORKDIR, reference_vcf,
        VARIANT_CALLERS=['strelka', 'mutect']):
    """Build variant density plots for all pipeline version and reference truth data."""
    output_dir = BENCHMARKING_CONFIG['final_output_dir']
    reference_vcf = allel.vcf_to_dataframe(reference_vcf)
    create_density_plot_for_all_chromosomes(
        reference_vcf, output_dir, 'reference.jpg')
    for pipeline, description in BENCHMARKING_CONFIG['neoantigen_discovery_pipelines'].items():
        for variant_caller in VARIANT_CALLERS:
            tested_vcf = os.path.join(
                WORKDIR, description['workdir'], 'output',
                f'{variant_caller}.vcf'
            )
            tested_vcf = allel.vcf_to_dataframe(tested_vcf)
            create_density_plot_for_all_chromosomes(
                tested_vcf, output_dir, f'{pipeline}_{variant_caller}.jpg')