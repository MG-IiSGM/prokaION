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
BLUE =  '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'


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

    file_info = os.stat(file_name) # Retrieve the file into to check if has size > 0

    if not os.path.isfile(file_name) or file_info.st_size == 0:
        logger.info(RED + BOLD + 'File: %s not found or empty\n' % file_name + END_FORMATTING)
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
        if root == input_dir: # This only apply to parent folder, not subdirectories
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


def execute_subprocess(cmd, isShell = False, isInfo = False):
    """
    https://crashcourse.housegordon.org/python-subprocess.html
    https://docs.python.org/3/library/subprocess.html 
    Execute and handle errors with subprocess, outputting stderr instead of the subprocess CalledProcessError
    """

    logger.debug('')
    logger.debug(cmd)

    if cmd[0] == 'samtools' or cmd[0] == 'bwa' or cmd[0] == 'artic':
        prog = ' '.join(cmd[0:2])
        param = cmd [3:]
    else:
        prog = cmd[0]
        param = cmd[1:]

    try:
        command = subprocess.run(cmd, shell = isShell, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        if command.returncode == 0:
            logger.debug(GREEN + DIM + 'Program %s successfully executed' % prog + END_FORMATTING)
        else:
            logger.info(RED + BOLD + 'Command %s FAILED\n' % prog + END_FORMATTING + BOLD + 'with parameters: ' + END_FORMATTING + ' '.join(param) + '\n' + BOLD + 'EXIT-CODE: %d\n' % command.returncode + 'ERROR:\n' + END_FORMATTING + command.stderr.decode().strip())

        if isInfo:
            logger.info(command.stdout.decode().strip())
        else:
            logger.debug(command.stdout.decode().strip())

        logger.debug(command.stderr.decode().strip())

    except OSError as e:
        sys.exit(RED + BOLD + "Failed to execute program '%s': %s" % (prog, str(e)) + END_FORMATTING)


def check_reanalysis(output_dir):
    output_dir = os.path.abspath(output_dir)
    # group = output_dir.split('/')[-1]

    bam_dir = os.path.join(output_dir, 'Bam')
    vcf_dir = os.path.join(output_dir, 'VCF')
    vcfr_dir = os.path.join(output_dir, 'VCF_recal')
    gvcf_dir = os.path.join(output_dir, 'GVCF')
    gvcfr_dir = os.path.join(output_dir, 'GVCF_recal')
    cov_dir = os.path.join(output_dir, 'Coverage')
    table_dir = os.path.join(output_dir, 'Table')

    previous_files = [bam_dir, vcf_dir, gvcf_dir, gvcfr_dir]

    # Check how many folder exist
    file_exist = sum([os.path.exists(x) for x in previous_files]) # True = 1, False = 0

    # Handle reanalysis: First time; reanalysis or reanalysis with aditional samples
    if file_exist > 0: # Already analysed

        samples_analyzed = os.listdir(bam_dir)
        samples_analyzed = len([x for x in samples_analyzed if '.bai' not in x and 'bqsr' in x])

        samples_fastq = os.listdir(output_dir)
        samples_fastq = len([x for x in samples_fastq if x.endswith('fastq.gz')])

        if samples_analyzed >= samples_fastq:
            logger.info(MAGENTA + '\nPrevious analysis detected, no new sequences added\n' + END_FORMATTING)

        else:
            logger.info(MAGENTA + '\nPrevious analysis detected, new sequences added\n' + END_FORMATTING)
            for root, _, files in os.walk(output_dir):
                if root == gvcf_dir or root == gvcfr_dir or root == vcfr_dir:
                    for name in files:
                        filename = os.path.join(root, name)
                        if (('GVCF_recal' in filename) or ('/VCF_recal' in filename)) and 'cohort' in filename and samples_analyzed < 100:
                            os.remove(filename)
                        elif 'cohort' in filename and '/GVCF/' in filename:
                            os.remove(filename)
                elif root == vcf_dir or root == table_dir:
                    for name in files:
                        filename = os.path.join(root, name)
                        if 'cohort' in filename or filename.endswith('.bed') or filename.endswith('.tab'):
                            os.remove(filename)
                elif root == cov_dir:
                    for name in files:
                        filename = os.path.join(root, name)
                        if 'coverage.tab' in filename:
                            os.remove(filename)
                        if 'poorly_covered.bed' in filename and samples_analyzed < 100:
                            os.remove(filename)


