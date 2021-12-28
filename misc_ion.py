#!/usr/bin/env python

import os
import sys
import re
import subprocess
import shutil
import logging
import datetime
import pandas as pd
import numpy as np
from statistics import mean
import multiprocessing


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

# Colors and formatting

"""
http://ozzmaker.com/add-colour-to-text-in-python/
The above ANSI escape code will set the text colour to bright green. The format is;
\033[  Escape code, this is always the same
1 = Style, 1 for normal.
32 = Text colour, 32 for bright green.
40m = Background colour, 40 is for black.
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


### Executing functions ###

def execute_subprocess(cmd, isShell=False, isInfo=False):
    """
    https://crashcourse.housegordon.org/python-subprocess.html
    https://docs.python.org/3/library/subprocess.html 
    Execute and handle errors with subprocess, outputting stderr instead of the subprocess CalledProcessError
    """

    logger.debug('')
    logger.debug(cmd)

    if cmd[0] == 'samtools' or cmd[0] == 'bwa' or cmd[0] == 'artic':
        prog = ' '.join(cmd[0:2])
        param = cmd[3:]
    else:
        prog = cmd[0]
        param = cmd[1:]

    try:
        command = subprocess.run(
            cmd, shell=isShell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if command.returncode == 0:
            logger.debug(
                GREEN + DIM + 'Program %s successfully executed' % prog + END_FORMATTING)
        else:
            logger.info(RED + BOLD + 'Command %s FAILED\n' % prog + END_FORMATTING + BOLD + 'with parameters: ' + END_FORMATTING + ' '.join(
                param) + '\n' + BOLD + 'EXIT-CODE: %d\n' % command.returncode + 'ERROR:\n' + END_FORMATTING + command.stderr.decode().strip())

        if isInfo:
            logger.info(command.stdout.decode().strip())
        else:
            logger.debug(command.stdout.decode().strip())

        logger.debug(command.stderr.decode().strip())

    except OSError as e:
        sys.exit(RED + BOLD + "Failed to execute program '%s': %s" % (prog,
                 str(e)) + END_FORMATTING)


### Manipulation of files and paths ###

def check_create_dir(path):
    # exists = os.path.isfile(path)
    # exists = os.path.isdir(path)

    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)


def check_file_exists(file_name):
    """
    Check file exist and is not 0Kb, if not program exit.
    """

    # Retrieve the file into to check if has size > 0
    file_info = os.stat(file_name)

    if not os.path.isfile(file_name) or file_info.st_size == 0:
        logger.info(RED + BOLD + 'File: %s not found or empty\n' %
                    file_name + END_FORMATTING)
        sys.exit(1)
    return os.path.isfile(file_name)


def check_remove_file(file_name):
    """
    Check file exist and remove it.
    """

    if os.path.exists(file_name):
        os.remove(file_name)


def extract_read_list(input_dir):

    input_dir = os.path.abspath(input_dir)
    all_files = []

    for root, _, files in os.walk(input_dir):
        if root == input_dir:  # This only apply to parent folder, not subdirectories
            for name in files:
                filename = os.path.join(root, name)
                is_files = re.match(r".*\.f(ast)*[aq5](\.gz)*", name)
                if is_files:
                    all_files.append(filename)

    all_files = sorted(all_files)

    return all_files


def extract_sample_list(file):

    basename_file = os.path.basename(file)
    basename_file = basename_file.split('.')[0]

    return basename_file


def check_reanalysis(output_dir, samples_to_analyze):

    output_dir = os.path.abspath(output_dir)

    new_samples = []

    variant_dir = os.path.join(output_dir, "Variants")
    compare_dir = os.path.join(output_dir, "Compare")

    previous_files = [variant_dir, compare_dir]

    # Check how many folders exist
    file_exist = sum([os.path.exists(x)
                     for x in previous_files])  # True = 1, False = 0

    # Handle reanalysis: First time; reanalysis or realysis with aditional samples
    if file_exist > 0:  # Already analysed

        previous_samples_list = os.listdir(variant_dir)

        if len(samples_to_analyze) == len(previous_samples_list):
            logger.info(
                MAGENTA + "\nPrevious analysis detected, no new sequences added\n" + END_FORMATTING)
        else:
            new_samples = set(samples_to_analyze) - set(previous_samples_list)
            logger.info(MAGENTA + "\nPrevious analysis detected, " +
                        str(len(new_samples)) + " new sequences added\n" + END_FORMATTING)

    return list(new_samples)


def file_to_list(file_name):

    list_F = []

    file_name_abs = os.path.abspath(file_name)
    with open(file_name_abs, 'r') as f:
        for line in f:
            list_F.append(line.strip())

    return list_F


### Processing Reference files ###

def optimal_chunk_size(reference):
    """
    https://www.bioinformatics.nl/emboss-explorer/
        http://www.sacs.ucsf.edu/Documentation/emboss/infoseq.html
    """

    input_reference = os.path.abspath(reference)
    input_folder = os.path.dirname(reference)
    chunks_size = 'chunks_size'
    chunks_size_path = os.path.join(input_folder, 'chunks_size')

    cmd_infoseq = ['infoseq', '-auto', '-only',
                   '-length', '-noheading', '-odirectory', input_folder, '-outfile', chunks_size, input_reference]
    # print(cmd_infoseq)
    execute_subprocess(cmd_infoseq, isShell=False)

    df = pd.read_csv(chunks_size_path, header=None)
    size = df.loc[df.index[0]]
    return size


def samtools_faidx(reference):
    # samtools faidx reference.fa

    input_reference = os.path.abspath(reference)
    fai_reference = input_reference + ".fai"

    if os.path.isfile(fai_reference):
        logger.info(fai_reference + " already EXISTS")
    else:
        logger.info(GREEN + "Indexing " + fai_reference + END_FORMATTING)
        cmd_faidx = 'samtools', 'faidx', reference
        execute_subprocess(cmd_faidx, isShell=False)


def create_reference_chunks(reference):

    ref_size = optimal_chunk_size(reference)
    min_freebayes_chunk_size = 1000
    chunk_cpu = multiprocessing.cpu_count()-2
    num_chunks = round(
        max(min_freebayes_chunk_size, int(ref_size) / chunk_cpu))

    input_reference = os.path.abspath(reference)
    input_folder = os.path.dirname(reference)
    out_reference_file = os.path.join(
        input_folder, 'reference.' + str(num_chunks) + '.regions')
    fai_reference = input_reference + '.fai'

    if os.path.isfile(out_reference_file):
        logger.info(out_reference_file + " already EXISTS")
    else:
        logger.info(GREEN + "Creating " + out_reference_file + END_FORMATTING)
        cmd_chunks = "fasta_generate_regions.py {} {} > {}".format(
            fai_reference, num_chunks, out_reference_file)
        execute_subprocess(cmd_chunks, isShell=True)

    return out_reference_file


### BAM Variant ###

def extract_indels(input_vcf):

    input_vcf = os.path.abspath(input_vcf)
    vcf_dir = ('/').join(input_vcf.split('/')[0:-1])
    output_indel_vcf = os.path.join(vcf_dir, 'snps.indel.vcf')

    with open(output_indel_vcf, 'w+') as fout:
        with open(input_vcf, 'r') as f:
            for line in f:
                if "TYPE=ins" in line or "TYPE=del" in line:
                    fout.write(line)


def merge_vcf(snp_vcf, indel_vcf):

    snp_vcf = os.path.abspath(snp_vcf)
    indel_vcf = os.path.abspath(indel_vcf)

    vcf_dir = ('/').join(snp_vcf.split('/')[0:-1])
    output_complete_vcf = os.path.join(vcf_dir, 'snps.all.vcf')

    with open(output_complete_vcf, 'w+') as fout:
        with open(snp_vcf, 'r') as f1:
            for line in f1:
                fout.write(line)
        with open(indel_vcf, 'r') as f2:
            for line in f2:
                if not line.startswith("#"):
                    fout.write(line)


def create_bamstat(input_bam, output_file, threads=36):

    cmd_bamstat = "samtools flagstat --threads {} {} > {}".format(
        str(threads), input_bam, output_file)
    # print(cmd_bamstat)
    execute_subprocess(cmd_bamstat, isShell=True)


def create_coverage(input_bam, output_file):

    cmd_coverage = "samtools depth -aa {} > {}".format(input_bam, output_file)
    # print(cmd_coverage)
    execute_subprocess(cmd_coverage, isShell=True)


### VCF processing ###

def import_VCF42_freebayes_to_tsv(vcf_file, sep='\t'):
    vcf_file = os.path.abspath(vcf_file)
    tsv_file = (".").join(vcf_file.split(".")[:-1]) + ".tsv"

    headers = []
    extra_fields = ['TYPE', 'DP', 'RO', 'AO']
    with open(tsv_file, 'w+') as fout:
        with open(vcf_file, 'r') as f:
            next_line = f.readline().strip()
            while next_line.startswith("#"):
                next_line = f.readline().strip()
                if next_line.startswith('#CHROM'):
                    headers = next_line.split('\t')

        headers = headers[:7] + extra_fields + ['OLDVAR']
        fout.write(("\t").join(headers) + "\n")

        with open(vcf_file, 'r') as f:
            for line in f:
                extra_field_list = []
                # and not 'complex' in line and not 'mnp' in line
                if not line.startswith("#"):
                    line_split = line.split(sep)[:8]
                    info = line_split[-1].split(";")
                    for field in extra_fields:
                        extra_field_list.append(
                            [x.split("=")[-1] for x in info if field in x][0])
                    if 'OLDVAR' in line:
                        extra_field_list.append([x.split("=")[-1]
                                                 for x in info if 'OLDVAR' in x][0].split(',')[0])
                    output_line = ("\t").join(
                        line_split[:7] + extra_field_list)
                    fout.write(output_line + "\n")


def import_tsv_freebayes_to_df(tsv_file, sep='\t'):

    tsv_file = os.path.abspath(tsv_file)

    df = pd.read_csv(tsv_file, sep=sep)

    df.rename(columns={'#CHROM': 'REGION', 'RO': 'REF_DP',
                       'DP': 'TOTAL_DP', 'AO': 'ALT_DP', 'QUAL': 'ALT_QUAL'}, inplace=True)

    df['REF_FREQ'] = df['REF_DP']/df['TOTAL_DP']
    df['ALT_FREQ'] = df['ALT_DP']/df['TOTAL_DP']

    df = df.sort_values(by=['POS']).reset_index(drop=True)

    return df[['REGION', 'POS', 'ID', 'REF', 'ALT', 'ALT_QUAL', 'FILTER', 'TOTAL_DP', 'TYPE', 'REF_DP', 'ALT_DP', 'REF_FREQ', 'ALT_FREQ', 'OLDVAR']]


def vcf_to_ivar_tsv(input_vcf, output_tsv):

    input_tsv = (".").join(input_vcf.split(".")[:-1]) + ".tsv"

    import_VCF42_freebayes_to_tsv(input_vcf)

    df = import_tsv_freebayes_to_df(input_tsv)
    df.to_csv(output_tsv, sep="\t", index=False)
