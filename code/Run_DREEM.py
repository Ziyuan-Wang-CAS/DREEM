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

Wrapper script for the 2 parts of the DREEM pipeline:
Step 1: Bit Vector
Step 2: EM Clustering

The 2 steps are run sequentially.
Mapping is done prior to Bit Vector if needed.

This version of the pipeline is designed to be run on a local machine.
"""
import os
import argparse
import time


def Run_DREEM():
    """
    """
    start_time = time.time()

    map_cmd = 'python3 Mapping.py {} {} {} {} {} {} {} {} {}'
    map_cmd = map_cmd.format(sample_name, ref_name, paired,
                             CPUS, L, X, input_dir, output_dir, picard_path)
    bv_cmd = 'python3 BitVector.py {} {} {} {} {} {} {} {} {} {} {} {}'
    bv_cmd = bv_cmd.format(sample_name, ref_name, START, END, SUR_BASES,
                           qscore_file, QSCORE_CUTOFF, input_dir, output_dir,
                           paired, picard_path, fastq)
    cluster_cmd = 'python3 EM_Clustering.py {} {} {} {} {} {} {} {} {} {} ' + \
                  '{} {} {} {} {} {} {}'
    cluster_cmd = cluster_cmd.format(sample_name, ref_name, START, END,
                                     MIN_ITS, INFO_THRESH, CONV_CUTOFF,
                                     NUM_RUNS, MAX_K, CPUS, NORM_PERC_BASES,
                                     exc_AC, SIG_THRESH, struct, input_dir,
                                     output_dir, ctrl)

    # Check if FASTQ option was specified. If so, run mapping
    if fastq:
        os.system(map_cmd)

    os.system(bv_cmd)
    os.system(cluster_cmd)

    end_time = time.time()
    time_taken = round((end_time - start_time) / 60, 2)
    print('Total time taken: {} mins'.format(time_taken))


if __name__ == '__main__':
    # Inputs from user
    parser = argparse.ArgumentParser(description='Run DREEM')
    parser.add_argument('input_dir', help='Path to input directory')
    parser.add_argument('output_dir', help='Path to output directory')
    parser.add_argument('sample_name', help='Name of sample')
    parser.add_argument('ref_name', help='Name of reference genome')
    parser.add_argument('START', help='Start pos in ref genome (1-based)')
    parser.add_argument('END', help='End pos in ref genome (1-based)')
    parser.add_argument('--single', help='Single-end sequencing',
                        action='store_true')
    parser.add_argument('--struct', help='Sec structure prediction',
                        action='store_true')
    parser.add_argument('--fastq', help='No BAM file',
                        action='store_true')
    parser.add_argument('--ctrl', help='Control sample',
                        action='store_true')
    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    sample_name = args.sample_name
    ref_name = args.ref_name
    START = int(args.START)
    END = int(args.END)
    single = args.single
    struct = args.struct
    fastq = args.fastq
    ctrl = args.ctrl

    # ----- Inputs not specified by the user. Modify these as needed. ------- #

    # Inputs for Step 1 - Mapping
    picard_path = './picard.jar'  # Picard jar file in cur dir
    CPUS = 24  # Number of threads to use for alignment. Also used by clustering.
    L = 12  # Seed length for Bowtie2
    X = 1000  # Max fragment length for valid paired-end alignments

    # Inputs for Step 2 - Bit Vector creation
    qscore_file = './phred_ascii.txt'  # ASCII char - Q score map
    SUR_BASES = 10  # Bases surrounding a deletion on each side
    QSCORE_CUTOFF = 20  # Qscore cutoff for a valid base

    # Inputs for Step 3 - EM Clustering
    MIN_ITS = 300  # Min number of iterations per EM run
    INFO_THRESH = 1.0  # Threshold for informative bits
    CONV_CUTOFF = 0.5  # Diff in log like for convergence
    NUM_RUNS = 10  # Number of independent EM runs per K
    MAX_K = 3  # Max K to work on
    SIG_THRESH = 0.005  # Threshold to distinguish signal from noise
    NORM_PERC_BASES = 10  # Perc of bases to use for normalization
    exc_AC = True  # exclude As and Cs?

    # Make sure the 'input' directory (inside the 'project' directory)
    # contains the appropriate files inside them.
    # It is not necessary to manually create the 'output' directory.
    if not input_dir.endswith('/'):
        input_dir += '/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    paired = False if single else True

    Run_DREEM()
