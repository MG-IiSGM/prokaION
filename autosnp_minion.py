#!/usr/bin/env python

import os
import re
import logging
import argparse
import sys
import subprocess
import datetime
import gzip

# Local application imports

from misc_ion import check_create_dir, check_file_exists, extract_read_list, extract_sample_list, execute_subprocess


logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com)
Version = 0
Created: 22 November 2021

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

    parser = argparse.ArgumentParser(
        prog='autosnp_minion.py', description='Pipeline to Variant Calling from MinION sequencing')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input', dest='input_dir', metavar='input_directory',
                             type=str, required=True, help='REQUIRED. Input directory containing all fastq files')

    input_group.add_argument('-r', '--reference', metavar='reference',
                             type=str, required=True, help='REQUIRED. File to map against')

    input_group.add_argument('-s', '--sample', metavar='sample', type=str,
                             required=False, help='Sample to identify further files')

    input_group.add_argument('-L', '--sample_list', type=str, required=False,
                             help='Sample names to analyse only in the file supplied')

    input_group.add_argument('-p', '--primers', type=str,
                             required=False, help='Bed file including primers to trim')

    output_group = parser.add_argument_group(
        'Output', 'Required parameter to output results')

    output_group.add_argument('-o', '--output', type=str, required=True,
                              help='REQUIRED. Output directory to extract all results')

    output_group.add_argument('-C', '--noclean', required=False,
                              action='store_false', help='Clean unwanted files for standard execution')


def run_snippy(minion_fastq, reference, output_dir, sample, threads=30, minqual=10, minfrac=0.1, mincov=1):
    """
    https://github.com/tseemann/snippy
    USAGE
        snippy [options] --outdir <dir> --ref <ref> --R1 <R1.fq.gz> --R1 <R2.fq.gz>
        snippy [options] --outdir <dir> --ref <ref> --ctgs <contigs.fa>
        snippy [options] --outdir <dir> --ref <ref> --bam <reads.bam>
    """

    # --cpus: Maximum number of CPU cores to use
    # --outdir: Output folder
    # --prefix: Prefix for output files (default 'snps')
    # --minqual: Minumum QUALITY in VCF column 6
    # --mincov: Minimum site depth to for calling alleles
    # --minfrac: Minumum proportion for variant evidence
    # --ref: Reference genome. Supports FASTA, GenBank, EMBL (not GFF)
    # --se: Single-end reads

    prefix = os.path.join(output_dir, sample)

    cmd_snippy = ['snippy', '--cpus', str(threads), '--outdir', prefix, '--minqual', str(minqual), '--mincov', str(
        mincov), '--minfrac', str(minfrac), '--ref', reference, '--se', minion_fastq]

    print(cmd_snippy)
    execute_subprocess(cmd_snippy, isShell=False)


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
