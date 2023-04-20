"""Images of each pipeline to compare should have been built and pushed to AWS.
"""
from sklearn.metrics import classification_report, confusion_matrix
import wget
import os
import pandas as pd
import allel
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse
from tqdm import tqdm
from datetime import datetime
from pathlib import Path


import logging

def run(args):
    assess_performance(
        args.output_dir, args.reference, args.predictions,
        args.predictions_tag)
    build_variant_density_plots(
        args.output_dir, args.reference, args.predictions,
        args.predictions_tag)


def extract_variant_calls_compare_reference(
        tested_vcf, reference_vcf):
    """Extract variant calls coordinates from vcf files and compare results.
    Returns the results of the comparison. 0 if the variant was not called, 1 if"""
    tested_vcf_df = allel.vcf_to_dataframe(tested_vcf)
    tested_vcf_df=tested_vcf_df[tested_vcf_df['FILTER_PASS'] == True]
    reference_vcf_df = allel.vcf_to_dataframe(reference_vcf)
    tested_vcf_df['CHROM'] = tested_vcf_df['CHROM'].astype(str)
    reference_vcf_df['CHROM'] = reference_vcf_df['CHROM'].astype(str)
    tested_vcf_df['POS'] = tested_vcf_df['POS'].astype(int)
    reference_vcf_df['POS'] = reference_vcf_df['POS'].astype(int)
    merged_vcf_df = tested_vcf_df.merge(
        reference_vcf_df, on=['CHROM', 'POS'], how='outer',
        suffixes=('_tested', '_reference')
    )
    return compare(merged_vcf_df)

def compare(merged_vcf_df):
    """Assign presence/absence call for reference/prediction.
       presence = 1, absence = 0. 
        Example of data in a row:
            CHROM                     chr1
            POS                      14498
            ID_tested                    .
            REF_tested                   G
            ALT_1_tested                 A
            ALT_2_tested               NaN
            ALT_3_tested               NaN
            QUAL_tested                NaN
            FILTER_PASS_tested        True
            ID_reference               NaN
            REF_reference              NaN
            ALT_1_reference            NaN
            ALT_2_reference            NaN
            ALT_3_reference            NaN
            QUAL_reference             NaN
            FILTER_PASS_reference      NaN
    """
    variant_calls_pred = []
    variant_calls_true = []
    with tqdm(total=merged_vcf_df.shape[0]) as progress_bar:
        for variant_tup in tqdm(merged_vcf_df.iterrows()):
            variant = variant_tup[1]
           
            if variant['FILTER_PASS_tested'] and \
                type(variant['FILTER_PASS_tested']) == bool:
                variant_calls_pred.append(1)
            else:
                variant_calls_pred.append(0)
            if variant['FILTER_PASS_reference'] and \
                type(variant['FILTER_PASS_reference']) == bool:
                variant_calls_true.append(1)
            else:
                variant_calls_true.append(0)
            progress_bar.update(1)
    return variant_calls_pred, variant_calls_true


def build_performance_comparision_dataframe(
       reference, predictions, predictions_tag):
    """Build a dataframe with the performance of the different pipelines"""
    performance_list = []
    confusion_matrix_list = []

    for prediction in predictions:
        y_pred, y_true = extract_variant_calls_compare_reference(
            prediction, reference)
        target_names = ['not called', 'called']
        results_dict = classification_report(
            y_true, y_pred, target_names=target_names, output_dict=True)
        results_dict['pipeline'] = predictions_tag
        variant_caller = Path(prediction).stem
        results_dict['variant_caller'] = variant_caller

        tn, fp, fn, tp = confusion_matrix(
            y_true, y_pred).ravel()
        confusion_matrix_list.append(
            pd.DataFrame(
                {
                    'TN':[tn],
                    'FP':[fp],
                    'FN':[fn],
                    'TP':[tp],
                    'pipeline': predictions_tag, 'variant_caller':variant_caller
                }
            )
        )

        performance_list.append(pd.DataFrame(results_dict))
    return pd.concat(performance_list), pd.concat(confusion_matrix_list)


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
    """Create density plots for all chromosomes present in the vcf file."""
    for key in vcf.keys():
        if key.endswith('/variants/POS'):
            pos = allel.SortedIndex(vcf[key])
            plot_windowed_variant_density(
                pos, window_size=100000,
                filename=os.path.join(
                    output_dir, f'{key}-{filename}'
                ),
                title='Variant density'
            )


def build_variant_density_plots(output_dir, reference_vcf, predictions, tag
        ):
    """Build variant density plots for all pipeline version and reference truth data."""
    reference_vcf = allel.vcf_to_dataframe(reference_vcf)
    create_density_plot_for_all_chromosomes(
        reference_vcf, output_dir, 'reference.jpg')
    for prediction in predictions:
        variant_caller = os.path.basename(prediction)
        prediction_df = allel.vcf_to_dataframe(prediction)
        create_density_plot_for_all_chromosomes(
            prediction_df, output_dir, f'{tag}_{variant_caller}.jpg')


def assess_performance(output_dir, reference, predictions, predictions_tag):
    """Assess performance of the pipeline versions
    """
    performance_df, confusion_mat_df = build_performance_comparision_dataframe(
        reference, predictions, predictions_tag)
    date_time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    performance_df.to_excel(
        os.path.join(output_dir, f'performance_{date_time}.xlsx')
    )
    confusion_mat_df.to_excel(
        os.path.join(output_dir, f'confusion_matrix_{date_time}.xlsx'),
        index=False
    )

def main():
    logging.basicConfig(
        filename='performance.log', encoding='utf-8', level=logging.DEBUG,
        format='%(asctime)s:%(levelname)s - %(message)s'
    )
    logging.debug('Start Benchmarking debug')
    parser = argparse.ArgumentParser(
                    prog='Benchmak',
                    description='Assess performance between reference and prediction vcf files'
                )
    parser.add_argument('--reference', help='Reference vcf file')
    parser.add_argument(
        '--predictions', help='Predictions vcf files', nargs='+')
    parser.add_argument(
        '--predictions-tag', help='Tag to identify predictions')
    parser.add_argument('--output-dir', help='Output folder')

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()

