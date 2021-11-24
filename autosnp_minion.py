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

from misc_ion import check_create_dir, check_file_exists, extract_read_list, extract_sample_list, execute_subprocess, check_reanalysis, file_to_list


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

    input_group.add_argument('-t', '--threads', type=int, dest='threads', required=False,
                             default=30, help='Threads to use (30 threads by default)')

    output_group = parser.add_argument_group(
        'Output', 'Required parameter to output results')

    output_group.add_argument('-o', '--output', type=str, required=True,
                              help='REQUIRED. Output directory to extract all results')

    arguments = parser.parse_args()

    return arguments


# def run_snippy(input_sample_dir, reference, out_variant_dir, threads=30, minqual=10, minfrac=0.1, mincov=1):
#     """
#     https://github.com/tseemann/snippy
#     USAGE
#         snippy [options] --outdir <dir> --ref <ref> --R1 <R1.fq.gz> --R1 <R2.fq.gz>
#         snippy [options] --outdir <dir> --ref <ref> --ctgs <contigs.fa>
#         snippy [options] --outdir <dir> --ref <ref> --bam <reads.bam>
#     """

#     # --cpus: Maximum number of CPU cores to use
#     # --outdir: Output folder
#     # --prefix: Prefix for output files (default 'snps')
#     # --minqual: Minumum QUALITY in VCF column 6
#     # --mincov: Minimum site depth to for calling alleles
#     # --minfrac: Minumum proportion for variant evidence
#     # --ref: Reference genome. Supports FASTA, GenBank, EMBL (not GFF)
#     # --se: Single-end reads

#     for root, _, files in os.walk(input_sample_dir):
#         for name in files:
#             if 'HQ' in name:
#                 # print(name)
#                 minion_fastq = os.path.join(root, name)
#                 # print(HQ_filename)
#                 minion_out_variant = os.path.join(
#                     out_variant_dir, os.path.basename(minion_fastq.split('.')[0]))
#                 check_create_dir(minion_out_variant)

#                 cmd_snippy = ['snippy', '--cpus', str(threads), '--outdir', minion_out_variant, '--minqual', str(
#                     minqual), '--mincov', str(mincov), '--minfrac', str(minfrac), '--ref', reference, '--se', minion_fastq]

#                 print(cmd_snippy)
#                 execute_subprocess(cmd_snippy, isShell=False)


def minimap2_mapping(out_samples_filtered_dir, out_sorted_bam, reference, threads=30):
    """
    https://github.com/lh3/minimap2
        # Oxford Nanopore genomic reads
        minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam

    http://www.htslib.org/doc/samtools.html
    """

    # -a: Output in the SAM format
    # -x: Preset (always applied before other options; see minimap2.1 for details) []
    #    - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
    #    - map-hifi - PacBio HiFi reads vs reference mapping
    #    - ava-pb/ava-ont - PacBio/Nanopore read overlap
    #    - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    #    - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
    #    - sr - genomic short-read mapping
    # -t: Number of threads

    # -b: Output BAM
    # -S: Ignored (input format is auto-detected)
    # -F: Only include reads with none of the FLAGS in INT present
    # --threads: Number of additional threads to use

    for root, _, files in os.walk(out_samples_filtered_dir):
        for name in files:
            HQ_filename = os.path.join(root, name)
            filename_out = name.split('.')[0].split('_')[1]
            # print(filename_out)
            filename_bam_out = os.path.join(
                out_sorted_bam, filename_out + '.sort.bam')
            filename_bai_out = os.path.join(
                out_sorted_bam, filename_out + '.sort.bam.bai')
            # print(filename_bam_out)

            if os.path.isfile(filename_bai_out):
                logger.info(YELLOW + DIM + BOLD + filename_bam_out +
                            ' EXIST\nOmmiting filtering for sample ' + filename_out + '\n' + END_FORMATTING)
            else:
                logger.info(GREEN + 'Mapping sample ' +
                            filename_out + END_FORMATTING)
                cmd_minimap2 = 'minimap2 -ax map-ont {} {} | samtools view -bS -F 4 - | samtools sort -o {}'.format(
                    reference, HQ_filename, filename_bam_out)
                # print(cmd_minimap2)
                execute_subprocess(cmd_minimap2, isShell=True)

                cmd_indexing = 'samtools', 'index', filename_bam_out
                # print(cmd_indexing)
                execute_subprocess(cmd_indexing, isShell=False)


if __name__ == '__main__':

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)
    in_samples_filtered_dir = os.path.join(
        input_dir, 'Samples_Fastq/Filtered_Fastq')
    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split('/')[-1]
    check_create_dir(output_dir)
    reference = os.path.abspath(args.reference)

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

    logger.info(
        "\n" + BLUE + '############### START VARIANT CALLING ###############' + END_FORMATTING + "\n")
    logger.info(args)

    # Obtain all fastq files from folder

    fastq = extract_read_list(in_samples_filtered_dir)

    # Check how many files will be analysed

    sample_list = []

    for sample in fastq:
        sample = extract_sample_list(sample)
        sample_list.append(sample)

    logger.info("\n" + CYAN + '{} Samples will be analysed: {}'.format(
        len(sample_list), ",".join(sample_list)) + END_FORMATTING)

    # Check if there are samples to filter out

    sample_list_F = []
    if args.sample_list == None:
        logger.info("\n" + 'No samples to filter')
        for sample in fastq:
            sample = extract_sample_list(sample)
            sample_list_F.append(sample)
    else:
        logger.info("Samples will be filtered")
        sample_list_F = file_to_list(args.sample_list)

    # Declare folders created in pipeline and key files

    out_bam_dir = os.path.join(output_dir, "Bam")
    check_create_dir(out_bam_dir)

    ############### START PIPELINE ###############

    # Mapping with minimap2, sorting Bam and indexing it

    minimap2_mapping(in_samples_filtered_dir, out_bam_dir,
                     reference=args.reference, threads=args.threads)

    logger.info("\n" + MAGENTA + BOLD +
                "##### END OF ONT VARIANT CALLING PIPELINE #####" + "\n" + END_FORMATTING)
