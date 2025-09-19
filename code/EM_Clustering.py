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

Step 2 of the DREEM pipeline: EM Clustering of bit vectors.
"""
import argparse
import time
import os
import BitVector_Functions
import EM_Plots
import EM_CombineRuns
import Run_EMJobs
import EM_Files
import sys
sys.setrecursionlimit(10000)

def EM_Clustering():
    """
    """
    print('Starting EM clustering...')

    for ref in refs_seq:  # Each seq in the ref genome

        start_time = time.time()

        bvfile_basename = '{}_{}_{}_{}'.format(sample_name, ref, START, END)
        outplot_dir = outfiles_dir + bvfile_basename + '/'
        if not os.path.exists(outplot_dir):
            os.makedirs(outplot_dir)
        else:  # Folder exists
            if os.path.exists(outplot_dir + 'log.txt'):  # Log file exists
                print('EM Clustering already done for', bvfile_basename)
                return

        wind_size = int(END) - int(START)
        norm_bases = int((wind_size * NORM_PERC_BASES) / 100)

        # Read the bit vector file and do the filtering
        input_file = output_dir + '/BitVector_Files/' + bvfile_basename + \
            '_bitvectors.txt'
        X = EM_Files.Load_BitVectors(input_file, INFO_THRESH, SIG_THRESH,
                                     exc_AC, output_dir, ctrl)

        K = 1  # Number of clusters
        cur_BIC = float('inf')  # Initialize BIC
        BIC_failed = False  # While test is not passed
        while not BIC_failed and K <= MAX_K:
            print('Working on K =', K)

            RUNS = NUM_RUNS if K != 1 else 1  # Only 1 Run for K=1
            ITS = MIN_ITS if K != 1 else 10  # Only 10 iters for K=1

            for run in range(1, RUNS + 1):
                print('Run number:', run)
                Run_EMJobs.Run_EMJob(X, bvfile_basename, ITS, INFO_THRESH,
                                     CONV_CUTOFF, SIG_THRESH,
                                     outplot_dir, K, CPUS, run)

            # Processing of results from the EM runs
            EM_CombineRuns.Post_Process(bvfile_basename, K, RUNS,
                                        cur_BIC, norm_bases, struct,
                                        input_dir, outplot_dir)

            # Check BIC
            latest_BIC = EM_CombineRuns.Collect_BestBIC(bvfile_basename, K,
                                                        outplot_dir)
            if latest_BIC > cur_BIC:  # BIC test has failed
                BIC_failed = True
            cur_BIC = latest_BIC  # Update BIC

            K += 1  # Move on to next K

        end_time = time.time()
        time_taken = round((end_time - start_time) / 60, 2)
        print('Time taken:', time_taken, 'mins')

        # Write params to log file
        EM_Plots.Log_File(bvfile_basename, NUM_RUNS, MIN_ITS,
                          CONV_CUTOFF, INFO_THRESH, SIG_THRESH, exc_AC,
                          norm_bases, K - 2, time_taken, outplot_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EM Clustering')
    parser.add_argument('sample_name', help='Name of sample')
    parser.add_argument('ref_name', help='Name of reference genome')
    parser.add_argument('START', help='Start coord for bit vector')
    parser.add_argument('END', help='End coord for bit vector')
    parser.add_argument('MIN_ITS', help='Min iterations per EM run')
    parser.add_argument('INFO_THRESH', help='Threshold for informative bits')
    parser.add_argument('CONV_CUTOFF', help='Diff in log like for convergence')
    parser.add_argument('NUM_RUNS', help='Number of EM runs per K')
    parser.add_argument('MAX_K', help='Max K to work on')
    parser.add_argument('CPUS', help='Num processors to use')
    parser.add_argument('NORM_PERC_BASES', help='Perc of bases to use for norm')
    parser.add_argument('exc_AC', help='Exclude signal from As & Cs?')
    parser.add_argument('SIG_THRESH', help='Signal threshold')
    parser.add_argument('struct', help='Sec structure prediction?')
    parser.add_argument('input_dir', help='Directory with input files')
    parser.add_argument('output_dir', help='Directory with output files')
    parser.add_argument('ctrl', help='Control sample')
    args = parser.parse_args()
    sample_name = args.sample_name
    ref_name = args.ref_name
    START = int(args.START)
    END = int(args.END)
    MIN_ITS = int(args.MIN_ITS)
    INFO_THRESH = float(args.INFO_THRESH)
    CONV_CUTOFF = float(args.CONV_CUTOFF)
    NUM_RUNS = int(args.NUM_RUNS)
    MAX_K = int(args.MAX_K)
    CPUS = int(args.CPUS)
    NORM_PERC_BASES = int(args.NORM_PERC_BASES)
    exc_AC = args.exc_AC
    SIG_THRESH = float(args.SIG_THRESH)
    struct = args.struct
    input_dir = args.input_dir
    output_dir = args.output_dir
    ctrl = args.ctrl

    ref_file = input_dir + ref_name + '.fasta'
    refs_seq = BitVector_Functions.Parse_FastaFile(ref_file)  # Ref seqs

    outfiles_dir = output_dir + '/EM_Clustering/'
    if not os.path.exists(outfiles_dir):
        os.makedirs(outfiles_dir)

    EM_Clustering()
