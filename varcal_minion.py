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


# Local application imports

from misc_ion import check_create_dir, check_file_exists, extract_read_list, extract_sample_list, execute_subprocess, check_reanalysis, file_to_list, samtools_faidx, create_reference_chunks, extract_indels, merge_vcf, vcf_to_ivar_tsv, create_bamstat, create_coverage, obtain_group_cov_stats, obtain_overal_stats, ivar_consensus, replace_consensus_header


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

    input_group.add_argument('-i', '--input', dest='input_dir', metavar='Input_directory',
                             type=str, required=True, help='REQUIRED. Input directory containing all fastq files')

    input_group.add_argument('-s', '--sample', metavar='sample', type=str,
                             required=False, help='Sample to identify further files')

    input_group.add_argument('-L', '--sample_list', type=str, required=False,
                             help='Sample names to analyse only in the file supplied')

    input_group.add_argument('-p', '--primers', type=str,
                             required=False, help='Bed file including primers to trim')

    input_group.add_argument('-t', '--threads', type=int, dest='threads', required=False,
                             default=36, help='Threads to use (36 threads by default)')

    variant_group = parser.add_argument_group(
        'Variant Calling', 'Variant Calling parameters')

    variant_group.add_argument('-B', '--bayes', required=False, action='store_true',
                               help='Variant Calling is done with freebayes-parallel')

    variant_group.add_argument('-f', '--min_allele_frequency', type=int, dest='min_allele_frequency', required=False,
                               default=0.1, help='Minimum fraction of observations supporting an alternate allele. Default: 0.1')

    variant_group.add_argument('-q', '--min_base_quality', type=int, dest='min_quality', required=False,
                               default=10, help='Exclude alleles from analysis below threshold. Default: 10')

    variant_group.add_argument('-m', '--min_mapping_quality', type=int, dest='min_mapping', required=False,
                               default=60, help='Exclude alignments from analysis below threshold. Default: 60')

    variant_group.add_argument('-F', '--min_frequency', type=int, dest='min_frequency', required=False,
                               default=0.7, help='Minimum fraction of observations to call a base. Default: 0.7')

    variant_group.add_argument('-d', '--min_depth', type=int, dest='min_depth', required=False,
                               default=8, help='Minimum depth to call a base. Default: 8')

    reference_group = parser.add_argument_group(
        'Reference', 'Reference parameters')

    reference_group.add_argument('-r', '--reference', metavar='Reference',
                                 type=str, required=True, help='REQUIRED. File to map against')

    reference_group.add_argument('--ploidy', type=int, dest='ploidy', required=False,
                                 default=1, help='Sets the default ploidy for the analysis')

    output_group = parser.add_argument_group(
        'Output', 'Required parameter to output results')

    output_group.add_argument('-o', '--output', type=str, required=True,
                              help='REQUIRED. Output directory to extract all results')

    quality_group = parser.add_argument_group(
        'Quality parameters', 'Parameters for diferent Quality conditions')

    quality_group.add_argument('-c', '--coverage20', type=int, default=50, required=False,
                               help='Minimum percentage of coverage at 20X to clasify as uncovered. Default: 50%')

    quality_group.add_argument('-n', '--min_snp', type=int, required=False,
                               default=30, help='SNP number to pass quality threshold. Default: 30 HQ SNP')

    arguments = parser.parse_args()

    return arguments


# def run_snippy(input_sample_dir, reference, out_variant_dir, threads=36, minqual=10, minfrac=0.1, mincov=1):
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


def minimap2_mapping(HQ_filename, filename_bam_out, reference):
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

    cmd_minimap2 = 'minimap2 -ax map-ont {} {} | samtools view -bS -F 4 - | samtools sort -o {}'.format(
        reference, HQ_filename, filename_bam_out)
    # print(cmd_minimap2)
    execute_subprocess(cmd_minimap2, isShell=True)

    cmd_indexing = 'samtools', 'index', filename_bam_out
    # print(cmd_indexing)
    execute_subprocess(cmd_indexing, isShell=False)


def freebayes_variant(reference, filename_bam_out, output_vcf, threads=36, frequency=0.1, base_qual=7, map_qual=5):
    """
    https://github.com/freebayes/freebayes
        # Freebayes-parallel
        freebayes-parallel <(fasta_generate_regions.py {fai_reference} {chunks}) {threads} {args} > {output}
    """

    # region_file:
    # --haplotype_length: Allow haplotype calls with contiguous embedded matches of up to this length
    # --use-best-n-alleles: Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores
    # --min-alternate-count: Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position
    # --min-alternate-fraction: Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position
    # -p: Sets the default ploidy for the analysis to N
    # --min-coverage:
    # -q: Exclude alleles from analysis if their supporting base quality is less than Q
    # -m: Exclude alignments from analysis if they have a mapping quality less than Q
    # --strict-vcf: Generate strict VCF format (FORMAT/GQ will be an int)

    region_file = create_reference_chunks(
        reference)

    cmd_bayes = 'freebayes-parallel {} {} -f {} --haplotype-length 0 --use-best-n-alleles 1 --min-alternate-count 0 --min-alternate-fraction {} -p 1 --min-coverage 1 -q {} -m {} --strict-vcf {} > {}'.format(
        region_file, str(threads), reference, str(frequency), str(base_qual), str(map_qual), filename_bam_out, output_vcf)
    print(cmd_bayes)
    execute_subprocess(cmd_bayes, isShell=True)


def bcftool_filter(output_raw_vcf, output_vcf):
    """
    https://samtools.github.io/bcftools/bcftools.html
        # bcftools view: View, subset and filter VCF files by position and filtering expression
        # bcftools annotate: Add or remove annotations
    """

    # --include: include sites for which Expression is true
    # --remove: list of annotations to remove

    # cmd_bcf_view = "bcftools view --include 'QUAL >= 10 && FMT/DP >= 1 && (FMT/AO)/(FMT/DP) >= 0.1' {} -o {}".format(
    #     output_raw_vcf, output_vcf)
    # # print(cmd_bcf_view)
    # execute_subprocess(cmd_bcf_view, isShell=True)

    # cmd_bcf_annot = 'bcftools annotate --remove ^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL {} -o {}'.format(
    #     output_vcf, output_vcf_filt)
    # # print(cmd_bcf_annot)
    # execute_subprocess(cmd_bcf_annot, isShell=True)

    cmd_bcf = "bcftools view --include 'QUAL >= 10 && FMT/DP >= 1 && (FMT/AO)/(FMT/DP) >= 0.1' {} | bcftools annotate --remove ^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL -o {}".format(
        output_raw_vcf, output_vcf)
    # print(cmd_bcf)
    execute_subprocess(cmd_bcf, isShell=True)


def snippy_sub(output_vcf_filt, output_vcf_sub):
    """
    https://github.com/tseemann/snippy/tree/master/bin
        # snippy-vcf_extract_subs: Convert MNP,COMPLEX into SNP and ignore INS,DEL
    """

    cmd_snippy_subs = 'snippy-vcf_extract_subs {} > {}'.format(
        output_vcf_filt, output_vcf_sub)
    # print(cmd_snippy_subs)
    execute_subprocess(cmd_snippy_subs, isShell=True)


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
        '\n' + BLUE + '############### START VARIANT CALLING ###############' + END_FORMATTING + '\n')
    logger.info(args)

    # Obtain all fastq files from folder

    fastq = extract_read_list(in_samples_filtered_dir)

    # Check how many files will be analysed

    sample_list = []

    for sample in fastq:
        sample = extract_sample_list(sample)
        # sample = sample.split('_')[1]
        sample_list.append(sample)

    # logger.info('\n' + CYAN + '{} Samples will be analysed: {}'.format(
    #     len(sample_list), ', '.join(sample_list)) + END_FORMATTING)

    # Check if there are samples to filter out

    sample_list_F = []
    if args.sample_list == None:
        logger.info('\n' + 'No samples to filter' + '\n')
        for sample in fastq:
            sample = extract_sample_list(sample)
            sample_list_F.append(sample)
    else:
        logger.info('Samples will be filtered')
        sample_list_F = file_to_list(args.sample_list)

    new_samples = check_reanalysis(args.output, sample_list_F)

    logger.info(CYAN + '\n%d samples will be analysed: %s' %
                (len(sample_list_F), ','.join(sample_list_F)) + END_FORMATTING)
    logger.info(CYAN + '\n%d NEW samples will be analysed: %s' %
                (len(new_samples), ','.join(new_samples)) + END_FORMATTING)

    # Declare folders created in pipeline and key files

    out_bam_dir = os.path.join(output_dir, 'Bam')
    check_create_dir(out_bam_dir)

    out_variant_dir = os.path.join(output_dir, 'Variants')
    check_create_dir(out_variant_dir)

    out_consensus_dir = os.path.join(output_dir, 'Consensus')
    check_create_dir(out_consensus_dir)

    out_stats_dir = os.path.join(output_dir, 'Stats')
    check_create_dir(out_stats_dir)
    out_stats_bamstats_dir = os.path.join(
        out_stats_dir, 'Bamstats')  # subfolder
    check_create_dir(out_stats_bamstats_dir)
    out_stats_coverage_dir = os.path.join(
        out_stats_dir, 'Coverage')  # subfolder
    check_create_dir(out_stats_coverage_dir)

    out_compare_dir = os.path.join(output_dir, 'Compare')
    check_create_dir(out_compare_dir)

    out_annot_dir = os.path.join(output_dir, 'Annotation')
    check_create_dir(out_annot_dir)
    out_annot_snpeff_dir = os.path.join(out_annot_dir, 'snpeff')  # subfolder
    check_create_dir(out_annot_snpeff_dir)
    out_annot_user_dir = os.path.join(out_annot_dir, 'user')  # subfolder
    check_create_dir(out_annot_user_dir)
    out_annot_user_aa_dir = os.path.join(out_annot_dir, 'user_aa')  # subfolder
    check_create_dir(out_annot_user_aa_dir)
    out_annot_blast_dir = os.path.join(out_annot_dir, 'blast')  # subfolder
    check_create_dir(out_annot_blast_dir)

    samtools_faidx(args.reference)

    create_reference_chunks(args.reference)

    ############### START PIPELINE ###############

    new_sample_number = 0

    for sample in fastq:
        # Extract sample name
        sample = extract_sample_list(sample)
        args.sample = sample
        if sample in sample_list_F:
            # Variant sample dir
            sample_variant_dir = os.path.join(out_variant_dir, sample)

            sample_number = str(sample_list_F.index(sample) + 1)
            sample_total = str(len(sample_list_F))

            if sample in new_samples:
                new_sample_number = str(int(new_sample_number) + 1)
                new_sample_total = str(len(new_samples))
                logger.info('\n' + WHITE_BG + 'STARTING SAMPLE: ' + sample + ' (' + sample_number + '/' +
                            sample_total + ')' + ' (' + new_sample_number + '/' + new_sample_total + ')' + END_FORMATTING)
            else:
                logger.info('\n' + WHITE_BG + 'STARTING SAMPLE: ' + sample +
                            ' (' + sample_number + '/' + sample_total + ')' + END_FORMATTING)

            # output_final_vcf = os.path.join(
            #     sample_variant_dir, 'snps.all.ivar.tsv')

    ##### MAPPING #####

    # Mapping with minimap2, sorting Bam and indexing it (also can be made with bwa index & bwa mem -x ont2d)

            # if not os.path.isfile(output_final_vcf):
            HQ_filename = os.path.join(
                in_samples_filtered_dir, sample + '.fastq.gz')
            # print(HQ_filename)
            # filename_out = sample.split('.')[0].split('_')[1]
            filename_out = sample
            # print(filename_out)
            filename_bam_out = os.path.join(
                out_bam_dir, filename_out + '.sort.bam')
            filename_bai_out = os.path.join(
                out_bam_dir, filename_out + '.sort.bam.bai')
            # print(filename_bam_out)

            logger.info(
                '\n' + GREEN + BOLD + 'STARTING ANALYSIS FOR SAMPLE ' + filename_out + END_FORMATTING)

            prior = datetime.datetime.now()

            if os.path.isfile(filename_bai_out):
                logger.info(YELLOW + filename_bam_out + ' EXIST\nOmmiting mapping for ' +
                            filename_out + END_FORMATTING)
            else:
                logger.info(GREEN + 'Mapping sample ' +
                            filename_out + END_FORMATTING)
                minimap2_mapping(
                    HQ_filename, filename_bam_out, reference=args.reference)

            after = datetime.datetime.now()
            print(('Done with function minimap2_mapping in: %s' %
                  (after - prior) + '\n'))

    ##### VARIANT CALLING #####

    # Variant calling with freebayes-parallel (also can be made with nanopolish, we should use nanopolish index & nanopolish variants)

            if filename_bam_out.endswith('bam'):
                # print(filename_bam_out)
                sample_variant_dir = os.path.join(
                    out_variant_dir, filename_out)
                # print(sample_variant_dir)
                check_create_dir(sample_variant_dir)
                output_raw_vcf = os.path.join(
                    sample_variant_dir, 'snps.raw.vcf')
                # print(output_raw_vcf)

                prior = datetime.datetime.now()

                if os.path.isfile(output_raw_vcf):
                    logger.info(
                        YELLOW + output_raw_vcf + ' EXIST\nOmmiting Variant Calling for sample ' + filename_out + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + 'Starting Variant Calling for sample ' + filename_out + END_FORMATTING)
                    freebayes_variant(args.reference, filename_bam_out, output_raw_vcf, threads=args.threads,
                                      frequency=args.min_allele_frequency, base_qual=args.min_quality, map_qual=args.min_mapping)

                after = datetime.datetime.now()
                print(('Done with function freebayes_variant in: %s' %
                      (after - prior) + '\n'))

    # Filtering the raw variant calling by quality, depth and frequency with bcftools. Also extracting complex variations and MNP with snippy-vcf_extract_subs

                output_vcf = os.path.join(
                    sample_variant_dir, 'snps.vcf')
                # print(output_vcf)
                output_vcf_filt = os.path.join(
                    sample_variant_dir, 'snps.filt.vcf')
                # print(output_vcf_filt)
                output_vcf_sub = os.path.join(
                    sample_variant_dir, 'snps.subs.vcf')
                # print(output_vcf_sub)

                prior = datetime.datetime.now()

                if os.path.isfile(output_vcf_sub):
                    logger.info(
                        YELLOW + output_vcf_sub + ' EXIST\nOmmiting Variant Calling filter in ' + filename_out + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + 'Variant Calling filtering in sample ' + filename_out + END_FORMATTING)
                    bcftool_filter(output_raw_vcf, output_vcf)
                    snippy_sub(output_vcf, output_vcf_sub)

                after = datetime.datetime.now()
                print(('Done with function bcftools_filter & snippy_sub in: %s' %
                      (after - prior) + '\n'))

    # Variant format combination, extracting INDELs and combining with subs.vcf (SNPs, MNPs and complex)

                out_variant_indel_sample = os.path.join(
                    sample_variant_dir, 'snps.indel.vcf')
                # print(out_variant_indel_sample)
                out_variant_all_sample = os.path.join(
                    sample_variant_dir, 'snps.all.vcf')
                # print(out_variant_all_sample)
                chrom_filename = os.path.join(
                    sample_variant_dir, 'snps.all.chromosome.vcf')
                # print(chrom_filename)

                prior = datetime.datetime.now()

                if os.path.isfile(out_variant_indel_sample):
                    logger.info(YELLOW + out_variant_indel_sample +
                                ' EXIST\nOmitting INDEL filtering in ' + filename_out + END_FORMATTING)
                else:
                    logger.info(GREEN + 'Filtering INDELs in ' +
                                filename_out + END_FORMATTING)
                    extract_indels(output_vcf)

                if os.path.isfile(out_variant_all_sample):
                    logger.info(YELLOW + out_variant_all_sample +
                                'EXIST\nOmitting VCF combination for sample ' + filename_out + END_FORMATTING)
                else:
                    logger.info(GREEN + 'Combining VCF in ' +
                                filename_out + END_FORMATTING)
                    merge_vcf(output_vcf_sub, out_variant_indel_sample)

                after = datetime.datetime.now()
                print(('Done with function extract_indels & merge_vcf in: %s' %
                      (after - prior) + '\n'))

    # Variant format adaptation

                out_variant_tsv_file = os.path.join(
                    sample_variant_dir, 'snps.all.ivar.tsv')
                # print(out_variant_tsv_file)

                prior = datetime.datetime.now()

                if os.path.isfile(out_variant_tsv_file):
                    logger.info(YELLOW + out_variant_tsv_file +
                                ' EXIST\nOmmiting format adaptation for ' + filename_out + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + 'Adapting variants format in sample ' + filename_out + END_FORMATTING)
                    vcf_to_ivar_tsv(out_variant_all_sample,
                                    out_variant_tsv_file)

                after = datetime.datetime.now()
                print(('Done with function vcf_to_ivar_tsv in: %s' %
                      (after - prior) + '\n'))

    ##### CONSENSUS #####

    # Building consensus fasta file with samtools mpileup and ivar consensus

                out_consensus_file = os.path.join(
                    out_consensus_dir, filename_out + '.fa')
                # print(out_consensus_file)

                prior = datetime.datetime.now()

                if os.path.isfile(out_consensus_file):
                    logger.info(YELLOW + out_consensus_file +
                                ' EXIST\nOmmiting Consensus for ' + filename_out + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + 'Creating Consensus in sample ' + filename_out + END_FORMATTING)
                    ivar_consensus(filename_bam_out, out_consensus_dir, filename_out, min_quality=args.min_quality,
                                   min_frequency_threshold=args.min_frequency, min_depth=8, uncovered_character='N')
                    replace_consensus_header(out_consensus_file)

                after = datetime.datetime.now()
                print(('Done with function ivar_consensus & replace_consensus_header in: %s' %
                      (after - prior) + '\n'))

    ##### CREATE STATS AND QUALITY FILTERS #####

    # Create Bamstats

        out_bamstats_name = filename_out + '.bamstats'
        out_bamstats_file = os.path.join(
            out_stats_bamstats_dir, out_bamstats_name)
        # print(out_bamstats_file)
        # print(filename_bam_out)

        prior = datetime.datetime.now()

        if os.path.isfile(out_bamstats_file):
            logger.info(YELLOW + out_bamstats_file +
                        'EXIST\nOmmiting Bamstats for ' + filename_out + END_FORMATTING)
        else:
            logger.info(GREEN + 'Creating Bamstats in sample ' +
                        filename_out + END_FORMATTING)
            create_bamstat(
                filename_bam_out, out_bamstats_file, threads=args.threads)

        after = datetime.datetime.now()
        print(('Done with function create_bamstat in: %s' %
              (after - prior) + '\n'))

    # Create Coverage

        out_coverage_name = filename_out + '.cov'
        out_coverage_file = os.path.join(
            out_stats_coverage_dir, out_coverage_name)
        # print(out_coverage_file)

        prior = datetime.datetime.now()

        if os.path.isfile(out_coverage_file):
            logger.info(YELLOW + out_coverage_file +
                        ' EXIST\nOmmiting Coverage for ' + filename_out + END_FORMATTING)
        else:
            logger.info(GREEN + 'Creating Coverage in sample ' +
                        filename_out + END_FORMATTING)
            create_coverage(filename_bam_out, out_coverage_file)

        after = datetime.datetime.now()
        print(('Done with function create_coverage in: %s' %
              (after - prior) + '\n'))

    # Coverage Output summary

    prior = datetime.datetime.now()

    logger.info(GREEN + BOLD + 'Creating summary report for coverage results in group ' +
                group_name + END_FORMATTING)
    obtain_group_cov_stats(out_stats_dir, group_name)

    # Reads and Variants output summary

    logger.info(GREEN + BOLD + 'Creating overal summary report in group ' +
                group_name + END_FORMATTING)
    obtain_overal_stats(out_stats_dir, output_dir, group_name)

    after = datetime.datetime.now()
    print(('Done with function obtain_group_cov_stats & obtain_overal_stats in: %s' %
          (after - prior) + '\n'))

    # Remove Uncovered

    logger.info('\n' + MAGENTA + BOLD +
                '##### END OF ONT VARIANT CALLING PIPELINE #####' + '\n' + END_FORMATTING)