#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# The MIT License (MIT)
# Copyright (c) <2019> <The Whitehead Institute for Biomedical Research>

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
Created on Tue Jul 30 2019

@author: harish

Initial step of the DREEM pipeline: Align reads in a FASTQ file to a ref genome.
This is done before Bit Vector if needed.

Output: SAM file. This and other files, such as QC stats files,
are created in separate output subdirectories.

There are 4 steps in this pipeline:
Step 1: QC of FASTQ files using FASTQC
Step 2: Adapter and quality trimming using TrimGalore (& another FASTQC check)
Step 3: Mapping using Bowtie2
Step 4: Post-mapping QC using Picard

Dependencies: FASTQC, TrimGalore, Bowtie2 and Picard.
Make sure you have installed all these tools and added them to your
computer's PATH environment variable
"""
import os
import argparse
import time
import datetime


def Map():
    """
    Perform QC, trimming and mapping for a given sample
    """

    start_time = time.time()
    sample_input_path = input_dir + sample_name
    sample_outfiles_path = outfiles_dir + sample_name + '_' + ref_name
    sample_outplots_path = outplots_dir + sample_name + '_' + ref_name
    bam_filename = sample_outfiles_path + '.bam'  # Bowtie2 output

    # ----------------  Step 1: QC using FASTQC ----------------------------- #

    mate1 = sample_input_path + '_mate1.fastq' if paired else \
        sample_input_path + '.fastq'
    mate2 = sample_input_path + '_mate2.fastq' if paired else ''
    if not os.path.exists(mate1):
        print('FastQ file ' + mate1 + ' does not exist.')
        return
    if paired:
        if not os.path.exists(mate2):
            print('FastQ file ' + mate2 + ' does not exist.')
            return
    if paired:
        fastqc_command = 'fastqc --extract {} {} --outdir={}'
        fastqc_command = fastqc_command.format(mate1, mate2, outplots_dir)
    else:
        fastqc_command = 'fastqc --extract {} --outdir={}'
        fastqc_command = fastqc_command.format(mate1, outplots_dir)
    os.system(fastqc_command)

    # ---------------- Step 2: Trimming using Trimgalore -------------------- #

    # Assumption of ASCII+33 encoding of Phred quality scores
    if paired:
        trim_command = 'trim_galore --fastqc --paired {} {} -o {}'
        trim_command = trim_command.format(mate1, mate2, outplots_dir)
    else:
        trim_command = 'trim_galore --fastqc {} -o {}'
        trim_command = trim_command.format(mate1, outplots_dir)
    os.system(trim_command)

    # ----------------- Step 3: Mapping using Bowtie2 ----------------------- #

    if not os.path.exists(refgenome_fasta):
        print('Ref fasta file ' + refgenome_fasta + ' does not exist.')
        return
    if not os.path.exists(refgenome_basename + '.4.bt2'):
        print('Bowtie2 indexing of ref genome has not been done yet.')
        return

    # TrimGalore output file names - supplied to Bowtie2
    trimmed_mate1 = outplots_dir + sample_name + '_mate1_val_1.fq' if \
        paired else outplots_dir + sample_name + '_trimmed.fq'
    trimmed_mate2 = outplots_dir + sample_name + '_mate2_val_2.fq' if \
        paired else ''

    if paired:
        map_cmd = 'bowtie2 --local --no-unal --no-discordant ' + \
            '--no-mixed -X {} -L {} -p {} -x {} -1 {} -2 {} -S {}.sam'
        map_cmd = map_cmd.format(X, L, p, refgenome_basename,
                                 trimmed_mate1, trimmed_mate2,
                                 sample_outfiles_path)
    else:
        map_cmd = 'bowtie2 --local --no-unal -L {} -p {} -x {} -U {} ' + \
            '-S {}.sam'
        map_cmd = map_cmd.format(L, p, refgenome_basename,
                                 trimmed_mate1, sample_outfiles_path)
    os.system(map_cmd)

    # ------------- Step 4: Post-mapping QC using Picard -------------------- #

    sam_filename = sample_outfiles_path + '.sam'  # Bowtie2 output
    sortsam_filename = sample_outfiles_path + '_sorted.sam'  # Picard output

    # Convert the SAM file to a BAM file
    convert_command = 'java -jar {} SamFormatConverter I={} O={}'
    convert_command = convert_command.format(picard_path, sam_filename,
                                             bam_filename)
    os.system(convert_command)

    # Sort the sam file - for the sake of Picard
    sort_command = 'java -jar {} SortSam I={} O={} SORT_ORDER=coordinate'
    sort_command = sort_command.format(picard_path, sam_filename,
                                       sortsam_filename)
    os.system(sort_command)

    # Collect various alignment metrics
    metrics_cmd = 'java -jar {} CollectMultipleMetrics I={} O={} R={}'
    metrics_cmd = metrics_cmd.format(picard_path, sortsam_filename,
                                     sample_outplots_path, refgenome_fasta)
    os.system(metrics_cmd)

    # Collect quality yield metrics
    qualyield_cmd = 'java -jar {} CollectQualityYieldMetrics I={} ' + \
                    'O={}_qual_yield_metrics.txt'
    qualyield_cmd = qualyield_cmd.format(picard_path, sortsam_filename,
                                         sample_outplots_path)
    os.system(qualyield_cmd)

    # Delete the sorted SAM file
    os.system('rm ' + sortsam_filename)

    now = datetime.datetime.now()
    end_time = time.time()
    time_taken = round((end_time - start_time) / 60, 2)

    # Write input parameters to log file
    log_file_name = sample_outplots_path + '_log.txt'
    log_file = open(log_file_name, 'w')
    log_file.write('Sample: ' + sample_name + '\n')
    log_file.write('Reference genome: ' + refgenome_fasta + '\n')
    log_file.write('Time taken: ' + str(time_taken) + ' mins\n')
    log_file.write('Finished at: ' + now.strftime('%Y-%m-%d %H:%M') + '\n')
    log_file.close()

    print('Finished mapping.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Align a FASTQ file to a ' +
                                     'reference genome')
    parser.add_argument('sample_name', help='Name of sample')
    parser.add_argument('ref_name', help='Name of reference genome')
    parser.add_argument('paired', help='Paired-end sequencing?')
    parser.add_argument('p', help='Number of threads for alignment')
    parser.add_argument('L', help='Seed length for Bowtie2')
    parser.add_argument('X', help='Max frag length for valid paired-end align')
    parser.add_argument('input_dir', help='Directory with input files')
    parser.add_argument('output_dir', help='Directory with output files')
    parser.add_argument('picard_path', help='Path to Picard jar file')
    args = parser.parse_args()
    sample_name = args.sample_name
    ref_name = args.ref_name
    paired = args.paired
    p = args.p
    L = args.L
    X = args.X
    input_dir = args.input_dir
    output_dir = args.output_dir
    picard_path = args.picard_path

    paired = True if paired == 'True' else False

    # Specify reference genome index files base name and fasta file
    # Index files (in same folder as the fasta file) are created by running:
    # $bowtie2-build 'refgenome.fasta' 'refgenome'
    refgenome_basename = input_dir + ref_name
    refgenome_fasta = refgenome_basename + '.fasta'  # Add extension

    # Output directories
    outfiles_dir = output_dir + '/Mapping_Files/'
    outplots_dir = output_dir + '/Mapping_Plots/'
    if not os.path.exists(outfiles_dir):
        os.makedirs(outfiles_dir)
    if not os.path.exists(outplots_dir):
        os.makedirs(outplots_dir)

    Map()
