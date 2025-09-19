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
"""
import argparse
import os
import EM_Plots
import EM_ExpandFold
import EM_ScatterClusters


def Collect_BestBIC(sample_name, K, outfiles_dir):
    """
    """
    K_dir = outfiles_dir + 'K_' + str(K) + '/'
    loglikes_file_name = K_dir + 'log_likelihoods.txt'
    loglikes_file = open(loglikes_file_name)
    for line in loglikes_file:
        line_split = line.strip().split()
        run_info = line_split[0]
        if run_info[-4:] == 'best':  # BIC from best run
            BIC = float(line_split[2])
            return BIC


def Post_Process(sample_name, K, RUNS, cur_BIC, norm_bases,
                 struct, input_dir, outfiles_dir):
    """
    """
    largest_loglike, BICs, log_likes, best_run = float('-inf'), [], [], ''

    for run in range(1, RUNS + 1):
        run_dir = outfiles_dir + '/K_' + str(K) + '/' + \
                      'run_' + str(run) + '/'

        largest_loglikefilename = run_dir + 'Largest_LogLike.txt'
        largest_loglikefile = open(largest_loglikefilename)
        log_like = float(largest_loglikefile.readline())
        largest_loglikefile.close()

        BIC_filename = run_dir + 'BIC.txt'
        BIC_file = open(BIC_filename)
        BIC = float(BIC_file.readline())
        BIC_file.close()

        log_likes.append(log_like)
        BICs.append(BIC)

        if log_like > largest_loglike:
            largest_loglike = log_like
            best_run = run

        os.system('rm ' + run_dir + 'Largest_LogLike.txt')
        os.system('rm ' + run_dir + 'BIC.txt')

    # Write to log likelihoods file
    EM_Plots.LogLikes_File(sample_name, K, RUNS, log_likes, BICs,
                           best_run, outfiles_dir)

    # Rename directory of best run
    orig_dir = outfiles_dir + 'K_' + str(K) + '/' + \
        'run_' + str(best_run) + '/'
    new_dir = outfiles_dir + 'K_' + str(K) + '/' + \
        'run_' + str(best_run) + '-best/'
    os.system('mv ' + orig_dir + ' ' + new_dir)

    clustmu_file = new_dir + 'Clusters_Mu.txt'

    # Folding with RNAstructure
    if struct == 'True':
        # Num bases on each side of region for secondary structure prediction
        num_bases = [0, 50, 100, 150, 200]
        for num_base in num_bases:
            EM_ExpandFold.ConstraintFoldDraw(input_dir, clustmu_file,
                                             num_base, num_base, norm_bases)

    # Scatter plot of reactivities
    if K > 1:
        EM_ScatterClusters.Scatter_Clusters(input_dir, clustmu_file)
