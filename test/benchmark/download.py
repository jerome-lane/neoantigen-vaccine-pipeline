"""Images of each pipeline to compare should have been built and pushed to AWS.
"""
import wget
import os
import argparse
from config import BENCHMARKING_CONFIG, SRA_TOOL_KIT_URL

WORKDIR = BENCHMARKING_CONFIG['workdir']

import logging


def run():
    sample_info = BENCHMARKING_CONFIG['sample_info']
    url = sample_info['SEQC2']['reference']['url']
    files = sample_info['SEQC2']['reference']['files']

    WORKDIR = BENCHMARKING_CONFIG['workdir']
    output_directory = os.path.join(
        WORKDIR, BENCHMARKING_CONFIG['SEQC2_reference_variant_calls_dir'])
    os.makedirs(output_directory, exist_ok=True)
    for file in files:
        file_url = os.path.join(url, file)
        destination = os.path.join(output_directory, file)
        if os.path.exists(destination):
            logging.debug('File %s already exists' % destination)
        else:
            logging.debug('Download %s to %s' % (file_url, destination))
            wget.download(file_url, out=destination)
    
    basename = os.path.basename(SRA_TOOL_KIT_URL)
    destination = os.path.join(WORKDIR, basename)
    if os.path.exists(destination):
            logging.debug('File %s already exists' % destination)
    else:
        logging.debug('Download %s' % (destination))
        wget.download(SRA_TOOL_KIT_URL, out=destination)

        os.system(
            f'tar -vxzf {WORKDIR}/sratoolkit.3.0.2-ubuntu64.tar.gz -C {WORKDIR} && export PATH=$PATH:{WORKDIR}/sratoolkit.3.0.2-ubuntu64/bin'
        )

    os.chdir(WORKDIR)

    for sample in sample_info['input']['samples']:
        sample_id = sample['ID']
        sample_folder = os.path.join(WORKDIR, sample_id)
        if not os.path.exists(sample_folder):
            os.system(
                f'{WORKDIR}/sratoolkit.3.0.2-ubuntu64/bin/prefetch {sample_id} --max-size 420000000000'
                )
            os.system(
                f'{WORKDIR}/sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump {sample_id}'
                )
            os.system(
                f'gzip {sample_id}')

def main():
    logging.basicConfig(
        filename='benchmark-download.log', encoding='utf-8',
        level=logging.DEBUG,
        format='%(asctime)s:%(levelname)s - %(message)s'
    )
    logging.debug('Start Benchmarking debug')
    parser = argparse.ArgumentParser(
                    prog='Benchmak-Download',
                    description='Download test set for benchmarking different version of openvax neoantigen discovery pipeline'
                )
    parser.parse_args()
    run()


if __name__ == '__main__':
    main()
