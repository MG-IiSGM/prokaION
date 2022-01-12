#!/usr/bin/env python

import os
import re
import logging
import argparse
import sys
import subprocess
import datetime
import gzip
import multiprocessing
import pandas as pd
import numpy as np
from pandarallel import pandarallel
import concurrent.futures
import multiprocessing
from sklearn.metrics import pairwise_distances, accuracy_score
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import scipy.cluster.hierarchy as shc
import scipy.spatial.distance as ssd  # pdist


logger = logging.getLogger()


END_FORMATTING = "\033[0m"
WHITE_BG = "\033[0;30;47m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
RED = "\033[31m"
GREEN = "\033[32m"
MAGENTA = "\033[35m"
BLUE = "\033[34m"
CYAN = "\033[36m"
YELLOW = "\033[93m"
DIM = "\033[2m"


def get_arguments():

    parser = argparse.ArgumentParser(
        prog="compare_snp_prokaion.py",
        description="Pipeline to compare variants (SNVs) with any non model organism. Specialised in Mycobacterium tuberculosis",
    )

    arguments = parser.parse_args()

    return arguments


def import_tsv_variants(
    tsv_file, sample, min_total_depth=4, min_alt_dp=7, only_snp=True
):

    input_file = os.path.abspath(tsv_file)

    df = pd.read_csv(input_file, sep="\t")
    df = df.drop_duplicates(subset=["POS", "REF", "ALT"], keep="first")

    df = df[((df.TOTAL_DP >= min_total_depth) & (df.ALT_DP >= min_alt_dp))]
    df = df[["REGION", "POS", "REF", "ALT", "ALT_FREQ", "TYPE"]]
    df = df.rename(columns={"ALT_FREQ": sample})

    if only_snp == True:
        df = df[df.TYPE == "snp"]
        df = df.drop(["TYPE"], axis=1)
        return df
    else:
        df = df.drop(["TYPE"], axis=1)
        return df


def extract_lowfreq(tsv_file, sample, min_total_depth=4, min_alt_dp=4, min_freq_include=0.7, only_snp=True):

    input_file = os.path.abspath(tsv_file)

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep='first')

    df = df[(df.ALT_DP <= min_alt_dp) & (df.ALT_FREQ >= min_freq_include)]

    df = df[['REGION', 'POS', 'REF', 'ALT', 'ALT_FREQ', 'TYPE']]
    df['ALT_FREQ'] = '?'
    df = df.rename(columns={'ALT_FREQ': sample})

    if only_snp == True:
        df = df[df.TYPE == 'snp']
        df = df.drop(['TYPE'], axis=1)
        return df
    else:
        df = df.drop(['TYPE'], axis=1)
        return df


def extract_uncovered(cov_file):

    base_file = os.path.basename(cov_file)
    input_file = os.path.abspath(cov_file)
    sample = base_file.split('.')[0]

    df = pd.read_csv(input_file, sep='\t', header=None)
    df.columns = ['REGION', 'POS', sample]
    df = df[df[sample] == 0]
    df = df.replace(0, '!')
    return df


def ddbb_create_intermediate(
    variant_dir,
    coverage_dir,
    min_freq_discard=0.1,
    min_alt_dp=7,
    only_snp=False,
    samples=False,
):

    df = pd.DataFrame(columns=["REGION", "POS", "REF", "ALT"])

    # Merge all raw
    for root, _, files in os.walk(variant_dir):
        for name in files:
            if name == "snps.all.ivar.tsv":
                sample = root.split("/")[-1]
                if samples == False:
                    logger.debug("Adding: " + sample)
                    filename = os.path.join(root, name)
                    dfv = import_tsv_variants(
                        filename,
                        sample,
                        min_total_depth=4,
                        min_alt_dp=4,
                        only_snp=only_snp,
                    )
                    df = df.merge(dfv, how="outer")
                else:
                    if sample in samples:
                        logger.debug("Adding: " + sample)
                        filename = os.path.join(root, name)
                        dfv = import_tsv_variants(
                            filename,
                            sample,
                            min_total_depth=4,
                            min_alt_dp=4,
                            only_snp=only_snp,
                        )
                        if dfv.shape[0] > 0:
                            df = df.merge(dfv, how="outer")
                        else:
                            logger.debug(sample + " HAS NO VARIANTS")
                            sample.remove(sample)

    # Round frequencies
    df = df.round(2)

    # Remove <= 0.1 (parameter in function)
    # IF HANDLE HETEROZYGOUS change this 0 for X or 0.5
    def handle_lowfreq(x): return None if x <= min_freq_discard else x
    df.iloc[:, 4:] = df.iloc[:, 4:].applymap(handle_lowfreq)

    # Drop all NaN rows
    df['AllNaN'] = df.apply(lambda x: x[4:].isnull().values.all(), axis=1)
    df = df[df.AllNaN == False]
    df = df.drop(['AllNaN'], axis=1).reset_index(drop=True)

    # Include poorly covered
    for root, _, files in os.walk(variant_dir):
        for name in files:
            if name == 'snps.all.ivar.tsv':
                filename = os.path.join(root, name)
                sample = root.split('/')[-1]
                if samples == False:
                    logger.debug('Adding lowfreqs: ' + sample)
                    dfl = extract_lowfreq(
                        filename, sample, min_total_depth=4, min_alt_dp=min_alt_dp, only_snp=only_snp)
                    df[sample].update(df[['REGION', 'POS', 'REF', 'ALT']].merge(
                        dfl, on=['REGION', 'POS', 'REF', 'ALT'], how='left')[sample])
                else:
                    if sample in samples:
                        logger.debug('Adding lowfreqs: ' + sample)
                        dfl = extract_lowfreq(
                            filename, sample, min_total_depth=4, min_alt_dp=min_alt_dp, only_snp=only_snp)
                        df[sample].update(df[['REGION', 'POS', 'REF', 'ALT']].merge(
                            dfl, on=['REGION', 'POS', 'REF', 'ALT'], how='left')[sample])

    indel_positions = df[(df['REF'].str.len() > 1) | (
        df['ALT'].str.len() > 1)].POS.tolist()
    indel_len = df[(df['REF'].str.len() > 1) | (
        df['ALT'].str.len() > 1)].REF.tolist()
    indel_len = [len(x) for x in indel_len]

    indel_positions_final = []

    for position, del_len in zip(indel_positions, indel_len):
        indel_positions_final = indel_positions_final + \
            [x for x in range(position - (5 + del_len),
                              position + (5 + del_len))]

    # Include uncovered
    samples_coverage = df.columns.tolist()[4:]
    for root, _, files in os.walk(coverage_dir):
        for name in files:
            if name.endswith('.cov'):
                filename = os.path.join(root, name)
                sample = name.split('.')[0]
                if sample in df.columns[4:]:
                    samples_coverage.remove(sample)
                    logger.debug('Adding uncovered: ' + sample)
                    dfc = extract_uncovered(filename)
                    dfc = dfc[~dfc.POS.isin(indel_positions_final)]
                    df[sample].update(df[['REGION', 'POS']].merge(
                        dfc, on=['REGION', 'POS'], how='left')[sample])

    if len(samples_coverage) > 0:
        logger.info('WARNING: ' + (',').join(samples_coverage) +
                    ' coverage file not found')

    # Asign 0 to rest (Absent)
    df = df.fillna(0)

    # Determine N (will help in poorly covered determination)

    def extract_sample_count(row):
        count_list = [i not in ['!', 0, '0'] for i in row[4:]]
        samples = np.array(df.columns[4:])
        # sample[np.array(count_list)] # Filter array with True/False array

        return (sum(count_list), (',').join(samples[np.array(count_list)]))

    if 'N' in df.columns:
        df = df.drop(['N', 'Samples'], axis=1)
    if 'Position' in df.columns:
        df = df.drop('Position', axis=1)

    df[['N', 'Samples']] = df.apply(
        extract_sample_count, axis=1, result_type='expand')

    df = df[df.N > 0]

    df['Position'] = df.apply(lambda x: ('|').join(
        [x['REGION'], x['REF'], str(x['POS']), x['ALT']]), axis=1)

    df = df.drop(['REGION', 'REF', 'POS', 'ALT'], axis=1)

    df = df[['Position', 'N', 'Samples'] +
            [col for col in df.columns if col not in ['Position', 'N', 'Samples']]]

    return df


def recheck_variant_rawvcf_intermediate(row, positions, alt_snps, variant_dir, min_cov_low_freq=7):
    """
    CU458896.1	3068036	.	G	A	262.784	.	AB=0.8;ABP=14.7363;AC=1;AF=0.5;AN=2;AO=12;CIGAR=1X;DP=15;DPB=15;DPRA=0;EPP=9.52472;EPPR=3.73412;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=1;NUMALT=1;ODDS=0.121453;PAIRED=0.333333;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=410;QR=108;RO=3;RPL=7;RPP=3.73412;RPPR=9.52472;RPR=5;RUN=1;SAF=4;SAP=5.9056;SAR=8;SRF=1;SRP=3.73412;SRR=2;TYPE=snp	GT:DP:AD:RO:QR:AO:QA:GL	0/1:15:3,12:3:108:12:410:-32.709,0,-5.55972
    """

    sample = row.name
    row_values = row.values
    return_list = [0] * len(positions)

    # Identify vcf
    variant_sample_folder = os.path.join(variant_dir, sample)
    filename = os.path.join(variant_sample_folder, "snps.raw.vcf")

    # Open file and retrieve output

    checked_positions = positions[:]

    with open(filename, 'r') as f:
        for line in f:
            # Parse line/variant in raw VCF
            if line.startswith('#CHROM'):
                headers = line.split('\t')
                # vcf_sample = headers[-1].strip()
            elif not line.startswith('#'):
                line_split = line.split('\t')
                vcf_position = line_split[1]
                params = line_split[-2].split(':')
                value_params = line_split[-1].split(':')
                depth_index = params.index('DP')

                # vcf_reference = line_split[0]
                vcf_alt_base = line_split[4]
                vcf_depth = int(value_params[depth_index])
                # ref_depth_index = params.index('RO')
                alt_depth_index = params.index('AO')

                try:
                    vcf_alt_depth = int(value_params[alt_depth_index])
                except:
                    value_params = [
                        int(x) for x in value_params[alt_depth_index].split(',')]
                    vcf_alt_depth = max(value_params)
                    vcf_alt_base = vcf_alt_base.split(',')[-1]

                vcf_alt_freq = round(vcf_alt_depth/vcf_depth, 2)

                # Process variant information for recalibration
                if vcf_position in positions:
                    if vcf_position in checked_positions:
                        checked_positions.remove(vcf_position)
                    position_index = positions.index(vcf_position)

                    if str(row[position_index]) == '0':
                        # vcf_ref_base = line_split[3]
                        position = positions[position_index]
                        alt_snp = alt_snps[position_index]

                        logger.debug('REC==>ORI:SAMPLE: {}\nPOS: {},  ALT: {}, DP: {}, FREQ: {}\nPOS: {}, ALT: {}, DP: {}, FREQ: {}'.format(
                            sample, vcf_position, vcf_alt_base, vcf_depth, vcf_alt_freq, position, alt_snp, vcf_alt_depth, row[position_index]))

                        try:
                            logger.debug('VALUE:\n{}, AO: {}'.format(
                                (':').join(value_params), alt_depth_index))
                        except:
                            logger.debug('problem in VALUES POS: {}, sample: {}'.format(
                                vcf_position, sample))

                        if vcf_depth <= min_cov_low_freq and vcf_depth > 0:
                            logger.debug('Position: {} LOWDEPTH: {}, DP: {}'.format(
                                vcf_position, vcf_alt_freq, vcf_depth))
                            row[position_index] = '?'
                        else:
                            if vcf_depth >= min_cov_low_freq and alt_snp == vcf_alt_base and vcf_alt_freq >= 0.1:
                                logger.debug('Position: {} RECOVERED from 0 to: {} in sample {}'.format(
                                    vcf_position, vcf_alt_freq, sample))
                                row[position_index] = vcf_alt_freq

                            # Handle INDEL incosistency between alt positions
                            elif (len(alt_snp) > 1) and vcf_alt_freq >= 0.1:
                                logger.debug('Position INDEL: {} RECOVERED from 0 to: {} in sample {}'.format(
                                    vcf_position, vcf_alt_freq, sample))
                                row[position_index] = vcf_alt_freq
                            elif alt_snp in vcf_alt_base:
                                logger.debug('Position COMPLEX: {} RECOVERED from 0 in {} to: {} in {} sample {}'.format(
                                    vcf_position, alt_snp, vcf_alt_freq, vcf_alt_base, sample))
                                row[position_index] = vcf_alt_freq

                elif 'complex' in line:
                    vcf_position = int(vcf_position)
                    positions_complex = [x for x in range(
                        vcf_position - 4, vcf_position + 4)]
                    positions_complex = [str(x) for x in positions_complex]
                    intersection = set(positions_complex).intersection(
                        set(checked_positions))
                    intersection = list(intersection)

                    if len(intersection) > 0:
                        for i in intersection:
                            if i in checked_positions:
                                checked_positions.remove(i)
                            position_index = positions.index(i)
                            if str(row[position_index]) == '0' or str(row[position_index]) == '!':
                                if vcf_depth <= min_cov_low_freq and vcf_depth > 0:
                                    logger.debug('Position: {} LOWDEPTH PRECOMPLEX: {}, DP: {}, sample {}'.format(
                                        i, vcf_alt_freq, vcf_depth, sample))
                                    row[position_index] = '?'
                                else:
                                    logger.debug('INTERSECTION PRECOMPLEX: Position: {}, AF: {}, sample: {}'.format(
                                        i, vcf_alt_freq, sample))
                                    row[position_index] = vcf_alt_freq
    return row


def recalibrate_ddbb_vcf_intermediate(snp_matrix_ddbb_file, variant_dir, min_cov_low_freq=7):
    """
    https://github.com/nalepae/pandarallel
    """

    pandarallel.initialize()

    selected_columns = ['Position', 'N',
                        'Samples', 'CHROM', 'REF', 'POS', 'ALT']

    df = pd.read_csv(snp_matrix_ddbb_file, sep='\t')
    df[['CHROM', 'REF', 'POS', 'ALT']
       ] = df.Position.str.split('|', expand=True)

    df = df[selected_columns +
            [col for col in df.columns if col not in selected_columns]]

    dft = df.transpose()

    samples = [str(x) for x in dft.index if x not in selected_columns]
    POSs = dft.loc['POS', :].tolist()
    # REFs = dft.loc['REF', :].tolist()
    ALTs = dft.loc['ALT', :].tolist()

    # n_samples = len(samples)

    header_df = dft.loc[selected_columns, :]
    samples_df = dft.drop(selected_columns, axis=0)

    samples_df = samples_df.parallel_apply(lambda x: recheck_variant_rawvcf_intermediate(
        x, POSs, ALTs, variant_dir, min_cov_low_freq=min_cov_low_freq), axis=1)

    final_dft = header_df.append(samples_df, ignore_index=True)

    final_df = final_dft.transpose()
    final_df.columns = selected_columns + samples
    final_df = final_df.drop(['CHROM', 'REF', 'POS', 'ALT'], axis=1)


def remove_position_range(df):

    INs = df[df['Position'].str.contains(
        r'\|[ATCG]{3,}\|[0-9]{1,}', regex=True)]

    bed_df = pd.DataFrame(columns=['#CHROM', 'REF', 'start', 'ALT'])
    bed_df[['#CHROM', 'REF', 'start', 'ALT']
           ] = INs['Position'].str.split('|', expand=True)
    bed_df['start'] = bed_df['start'].astype(int)

    # Include only positions in +1 position
    bed_df['start'] = bed_df['start'] + 1
    bed_df['lenREF'] = bed_df['REF'].str.len()

    # Commented to avoid collision with insertions
    bed_df['lenALT'] = bed_df['ALT'].str.len()
    bed_df['end'] = bed_df['start'] + bed_df['lenREF'] - 1

    for _, row in df.iterrows():
        position_number = int(row.Position.split('|')[2])
        if any(start <= position_number <= end for (start, end) in zip(bed_df.start.values.tolist(), bed_df.end.values.tolist())) and not re.search('\|[ATCG]{2,}', row.Position):
            logger.debug(
                'Position: {} removed found in INDEL'.format(row.Position))
            df = df[df.Position != row.Position]

    return df


def extract_complex_list(variant_dir, samples=False):

    all_complex = []

    for root, _, files in os.walk(variant_dir):
        for name in files:
            if name == 'snps.all.ivar.tsv':
                sample = root.split('/')[-1]
                filename = os.path.join(root, name)
                if samples == False:
                    logger.debug('Obtaining complex in: ' + sample)
                    df = pd.read_csv(filename, sep='\t')
                    sample_complex = df[~df.OLDVAR.isna()]['POS'].tolist()
                    sample_complex = [int(x) for x in sample_complex]
                    all_complex = all_complex + sample_complex
                else:
                    if sample in samples:
                        logger.debug('Obtaining complex in: ' + sample)
                        df = pd.read_csv(filename, sep='\t')
                        sample_complex = df[~df.OLDVAR.isna()]['POS'].tolist()
                        sample_complex = [int(x) for x in sample_complex]
                        all_complex = all_complex + sample_complex

    return sorted(set(all_complex))


def add_window_distance(vcf_df, window_size=10):
    """
    Add a column indicating the maximum number of SNPs in a windows of 10 or supplied distance
    """

    list_pos = vcf_df.POS.to_list()  # All positions
    set_pos = set(list_pos)  # To set for later comparing


def revised_df(df, out_dir=False, complex_pos=False, min_freq_include=0.7, min_threshold_discard_uncov_sample=0.6, min_threshold_discard_uncov_pos=0.6, min_threshold_discard_htz_sample=0.6, min_threshold_discard_htz_pos=0.6, min_threshold_discard_all_pos=0.6, min_threshold_discard_all_sample=0.6, remove_faulty=True, drop_samples=True, drop_positions=True, windows_size_discard=2):

    if remove_faulty == True:

        uncovered_positions = df.iloc[:, 3:].apply(lambda x:  sum(
            [i in ['!', '?'] for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=1)
        heterozygous_positions = df.iloc[:, 3:].apply(lambda x: sum([(i not in ['!', '?', 0, 1, '0', '1']) and (float(
            i) < min_freq_include) and (float(i) > 0.1) for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=1)

        report_position = pd.DataFrame({'Position': df.Position, 'uncov_fract': uncovered_positions,
                                       'htz_frac': heterozygous_positions, 'faulty_frac': uncovered_positions + heterozygous_positions})
        faulty_positions = report_position['Position'][(report_position.uncov_fract >= min_threshold_discard_uncov_pos) | (
            report_position.htz_frac >= min_threshold_discard_htz_pos) | (report_position.faulty_frac >= min_threshold_discard_all_pos)].tolist()

        # Calculate close SNPs/INDELs and remove those with 2 or more mutations in 10bp
        df['POS'] = df.apply(lambda x: x.Position.split('|')[2], axis=1)
        df['POS'] = df['POS'].astype(int)
        df = df.sort_values('POS')
        add_window_distance(df)
