"""Images of each pipeline to compare should have been built and pushed to AWS.
"""
from plots import build_variant_density_plots
from manage_reference_set import create_reference_vcf_file
from download import run as download_reference_files
from assess_performance import build_performance_comparison_dataframe
import os
import argparse
import pandas as pd
import glob
from datetime import datetime


import logging

from config import BENCHMARKING_CONFIG
WORKDIR = BENCHMARKING_CONFIG['workdir']
REFERENCE_VCF = os.path.join(
    WORKDIR, 
    BENCHMARKING_CONFIG['SEQC2_reference_variant_calls_dir'],
    'reference.vcf'
)
OLD_VARIANT_CALLER_IMAGE_NAME = BENCHMARKING_CONFIG['neoantigen_discovery_pipelines']['old']['image']
NEW_VARIANT_CALLER_IMAGE_NAME = BENCHMARKING_CONFIG['neoantigen_discovery_pipelines']['new']['image']
VARIANT_CALL_VCF_FILES_FOLDER = os.path.join(
    WORKDIR, BENCHMARKING_CONFIG['variant_call_vcf_copies']['workdir'])
FINAL_OUTPUT_DIR = os.path.join(
    WORKDIR, BENCHMARKING_CONFIG['final_output_dir'])

def run(args):
    sample_info = BENCHMARKING_CONFIG['sample_info']
    url = sample_info['SEQC2']['reference']['url']
    files = sample_info['SEQC2']['reference']['files']

    if args.download:
        download_reference_files(url, files, args.output_dir)
    if args.concatenate_reference_files:
        create_reference_vcf_file(WORKDIR, BENCHMARKING_CONFIG)
    if args.build_docker_images:
        os.system(
           f'docker build docker/Dockerfile -t {OLD_VARIANT_CALLER_IMAGE_NAME} .'
        )
        os.system(
           f'docker build docker/Dockerfile -t {NEW_VARIANT_CALLER_IMAGE_NAME} .'
        )
    if args.run_pipeline_to_benchmark:
        os.system('sh run_neoantigen_pipeline.sh')
    if args.assess_performance:
        perf_dfs = []
        conf_dfs = []
        for tag in ['new', 'old']:
            output_dir = os.path.join(VARIANT_CALL_VCF_FILES_FOLDER, tag)
            vcfs = glob.glob(f'{output_dir}/*.vcf')
            logging.debug(f'vcfs: {vcfs}')
            perf_df, conf_df = build_performance_comparison_dataframe(
                REFERENCE_VCF, vcfs, tag)
            perf_dfs.append(perf_df)
            conf_dfs.append(conf_df)
        
        date_time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        pd.concat(perf_dfs).to_csv(
            os.path.join(FINAL_OUTPUT_DIR, f'performance-all-{date_time}.csv')
        )
        pd.concat(conf_dfs).to_csv(
            os.path.join(FINAL_OUTPUT_DIR, f'confusion-matrix-all-{date_time}.csv')
        )
    if args.build_variant_density_plots:
        build_variant_density_plots(WORKDIR, REFERENCE_VCF)


def main():
    logging.basicConfig(
        filename='benchmarking.log', encoding='utf-8', level=logging.DEBUG,
        format='%(asctime)s:%(levelname)s - %(message)s'
    )
    logging.debug('Start Benchmarking debug')
    parser = argparse.ArgumentParser(
                    prog='Benchmak',
                    description='Benchmark different version of openvax neoantigen discovery pipeline'
                )
    parser.add_argument('-d', '--download', action='store_true')
    parser.add_argument('-r', '--run_pipeline_to_benchmark', action='store_true')

    parser.add_argument('-a', '--assess_performance', action='store_true')
    parser.add_argument('-c', '--concatenate_reference_files', action='store_true')
    parser.add_argument('-i', '--build_docker_images', action='store_true')

    parser.add_argument('-p', '--build_variant_density_plots', action='store_true')

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()

