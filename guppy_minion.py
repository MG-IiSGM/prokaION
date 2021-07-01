#!/usr/bin/env python

import os
import re
import logging
import argparse
import sys
import subprocess
import datetime

# Local application imports

from misc_ion import check_create_dir, check_file_exists, extract_read_list, extract_sample_list, execute_subprocess
from preprocess_ion import fastqc_quality


logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com), Pedro J. Sola (pedroscampoy@gmail.com)
Version = 0
Created: 24 March 2021

TODO:
    Adapt check_reanalysis
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE = '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def get_arguments():


    parser = argparse.ArgumentParser(prog = 'guppy_minion.py', description = 'Pipeline to basecalling and barcoding fast5 files from minION sequencing')

    parser.add_argument('-i', '--input', dest = 'input_dir', metavar = 'input_directory', type = str, required = True, help = 'REQUIRED. Input directory containing all fast5 files')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'REQUIRED. Output directory to extract all results')

    parser.add_argument('-s', '--samples', metavar = 'Samples', type = str, required = False, help = 'Sample list for conversion from barcode to samples ID')

    parser.add_argument('-c', '--config', type = str, default = 'dna_r9.4.1_450bps_fast.cfg', required = False, help = 'REQUIRED. Config parameter for guppy_basecalling. High-accuracy mode basecalling by default')

    parser.add_argument('-b', '--require_barcodes_both_ends', required = False, action = 'store_true', help = 'Require barcodes at both ends. By default it only requires the barcode at one end for the sequences identification')

    parser.add_argument('--barcode_kit', type = str, required = False, default = 'EXP-NBD196', help = 'Kit of barcodes used')

    parser.add_argument('-t', '--threads', type = int, dest = 'threads', required = False, default = 30, help = 'Threads to use (30 threads by default)')

    parser.add_argument('--num_callers', type = int, dest = 'num_callers', required = False, default = 10, help = 'Number of parallel basecallers')


    arguments = parser.parse_args()

    return arguments


def basecalling_ion(input_dir, out_basecalling_dir, config = 'dna_r9.4.1_450bps_fast.cfg', callers = 3, chunks = 2048, threads = 10):

    # -i: Path to input fast5 files
    # -s: Path to save fastq files
    # -c: Config file to use > https://community.nanoporetech.com/posts/guppy-v5-0-7-release-note (fast // hac // sup)
    # --num_callers: Number of parallel basecallers to Basecaller, if supplied will form part
    # --cpu_threads_per_caller: Number of CPU worker threads per basecaller
    # --chunks_per_runner: Maximum chunks per runner
    # --compress_fastq: Compress fastq output files with gzip

    cmd = ['guppy_basecaller', '-i', input_dir, '-s', out_basecalling_dir, '-c', config, '--num_callers', str(callers), '--chunks_per_runner', str(chunks), '--cpu_threads_per_caller', str(threads), '--compress_fastq']

    print(cmd)
    execute_subprocess(cmd, isShell = False)


def barcoding_ion(out_basecalling_dir, out_barcoding_dir, require_barcodes_both_ends = False, barcode_kit = 'EXP-NBD196', threads = 30):

    # -i: Path to input files
    # -r: Search for input file recursively
    # -s: Path to save files
    # -t: Number of worker threads
    # --fastq_out: Output Fastq files
    # --compress_fastq: Compress fastq output files with gzip
    # --barcode_kits: Space separated list of barcoding kit(s) or expansion kit(s) to detect against. Must be in double quotes
    # --require_barcodes_both_ends: Reads will only be classified if there is a barcode above the min_score at both ends of the read

    if require_barcodes_both_ends:
        logger.info(GREEN + BOLD + 'Barcodes are being used at both ends')
        require_barcodes_both_ends = "--require_barcodes_both_ends"
    else:
        logger.info(YELLOW + DIM + BOLD + 'Barcodes are being used on at least 1 of the ends')
        require_barcodes_both_ends = ""

    cmd = ['guppy_barcoder', '-i', out_basecalling_dir, '-s', out_barcoding_dir, '-r', require_barcodes_both_ends, '--barcode_kit', barcode_kit, '-t', str(threads), '--fastq_out', '--compress_fastq']

    print(cmd)
    execute_subprocess(cmd, isShell = False)


def read_filtering(out_barcoding_dir, out_samples_dir, summary = False, min_length = 270, max_length = 550):
    
    # --directory:
    # --prefix:
    # --min-length:
    # --max-length:
    # --output:

    with open(summary, 'r') as f:
        for line in f:
            barcode = os.path.join(out_barcoding_dir, line.split('\t')[0].strip())
            sample = line.split('\t')[1].strip()
            output_samples = os.path.join(out_samples_dir, sample + '.fastq')
            # print(output_samples)

            cmd = ['artic', 'guppyplex', '--directory', barcode, '--prefix', sample, '--min-length', str(min_length), '--max-length', str(max_length), '--output', output_samples]

            # print(cmd)
            execute_subprocess(cmd, isShell = False, isInfo = True)

            cmd_compress = ['bgzip', output_samples, '--threads', str(args.threads)]
            # print(cmd_compress)
            execute_subprocess(cmd_compress, isShell = False)



if __name__ == '__main__':

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split('/')[-1]
    check_create_dir(output_dir)

    # Logging
    # Create log file with date and time

    right_now = str(datetime.datetime.now())
    right_now_full = '_'.join(right_now.split(' '))
    log_filename = group_name + '_' + right_now_full + '.log'
    log_folder = os.path.join(output_dir, 'Logs')
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info("\n" + BLUE + '############### START PROCESSING FAST5 FILES ###############' + END_FORMATTING + "\n")
    logger.info(args)


    # Obtain all fast5 files from folder

    fast5 = extract_read_list(args.input_dir)


    # Check how many files will be analysed

    sample_list = []

    for sample in fast5:
        sample = extract_sample_list(sample)
        sample_list.append(sample)

    logger.info("\n" + CYAN + "{} Samples will be analysed: {}".format(len(sample_list), ",".join(sample_list)) + END_FORMATTING)


    # Declare folders created in pipeline and key files

    out_basecalling_dir = os.path.join(output_dir, 'Basecalling')
    check_create_dir(out_basecalling_dir)
    out_barcoding_dir = os.path.join(output_dir, 'Barcoding')
    check_create_dir(out_barcoding_dir)
    out_samples_dir = os.path.join(output_dir, 'Samples_Fastq')
    check_create_dir(out_samples_dir)
    out_fastqc_dir = os.path.join(output_dir, 'Quality')
    check_create_dir(out_fastqc_dir)


    ############### Start pipeline ###############

    # Basecalling

    logger.info("\n" + GREEN + "STARTING BASECALLING" + END_FORMATTING + "\n")

    basecalling_ion(input_dir, out_basecalling_dir, config = args.config, callers = args.num_callers, chunks = 2048, threads = args.threads)


    # Barcoding

    logger.info("\n" + GREEN + "STARTING BARCODING/DEMULTIPLEX" + END_FORMATTING + "\n")

    barcoding_ion(out_basecalling_dir, out_barcoding_dir, barcode_kit = args.barcode_kit, threads = args.threads, require_barcodes_both_ends = args.require_barcodes_both_ends)


    # Read Filtering

    logger.info("\n" + GREEN + "SAMPLE FILTERING" + END_FORMATTING + "\n")

    read_filtering(out_barcoding_dir, out_samples_dir, summary = args.samples)


    # Quality Check

    logger.info("\n" + GREEN + "QUALITY CHECK IN RAW" + END_FORMATTING)

    fastq = extract_read_list(out_samples_dir)

    sample_fastq = []

    for sample in fastq:
        sample_base = os.path.basename(sample)
        sample_fastq.append(os.path.basename(sample))

        out_qc_raw_fastq = ('.').join(sample.split('/')[-1].split('.')[0:-2]) + '_fastqc.html'
        out_qc_raw_file_fastq = os.path.join(out_fastqc_dir, out_qc_raw_fastq)

        if os.path.isfile(out_qc_raw_file_fastq):
            logger.info(YELLOW + DIM + out_qc_raw_file_fastq + ' EXIST\nOmmiting QC for sample ' + sample_base + END_FORMATTING)
        else:
            logger.info(GREEN + 'Checking quality in sample ' + sample_base + END_FORMATTING)
            fastqc_quality(sample, out_fastqc_dir, args.threads)


    logger.info("\n" + MAGENTA + BOLD +
                "#####END OF PIPELINE#####" + END_FORMATTING + "\n")