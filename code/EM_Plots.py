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

Does all the EM Clustering plots
"""
import os
import plotly
import plotly.graph_objs as go
from plotly import tools
import datetime


def Run_Plots(sample_name, X, K, log_like_list, final_mu, final_obs_pi,
              final_real_pi, resps, BIC, outplots_dir, run):
    """
    """
    K_dir = outplots_dir + '/K_' + str(K) + '/'
    if not os.path.exists(K_dir):
        os.makedirs(K_dir)
    run_dir = K_dir + 'run_' + str(run) + '/'
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)

    indices = X.indices.split(',')
    start, end = int(indices[0]), int(indices[1])
    seq = X.seq

    # File 1 - List of log likelihoods
    outfile_name1 = run_dir + 'Log_Likelihoods.txt'
    outfile1 = open(outfile_name1, 'w')
    for log_like in log_like_list:
        str_loglike = str(round(log_like, 2))
        outfile1.write(str_loglike + '\n')
    outfile1.close()

    # File 2 - Largest log likelihood
    outfile_name2 = run_dir + 'Largest_LogLike.txt'
    outfile2 = open(outfile_name2, 'w')
    str_loglike = str(round(log_like_list[-1], 2))
    outfile2.write(str_loglike + '\n')
    outfile2.close()

    # File 3 - BIC
    outfile_name3 = run_dir + 'BIC.txt'
    outfile3 = open(outfile_name3, 'w')
    str_BIC = str(round(BIC, 2))
    outfile3.write(str_BIC + '\n')
    outfile3.close()

    # File 4 - Cluster mus
    outfile_name4 = run_dir + 'Clusters_Mu.txt'
    outfile4 = open(outfile_name4, 'w')
    outfile4.write('@ref' + '\t' + X.ref_file + ';' + X.ref + '\t' +
                   seq[start - 1:end] + '\n')
    outfile4.write('@coordinates:length' + '\t' + str(start) + ',' +
                   str(end) + ':' + str(end - start + 1) + '\n')
    outfile4.write('Position')
    for i in range(len(final_mu)):
        outfile4.write('\tCluster_' + str(i + 1))
    outfile4.write('\n')
    for i in range(start, end + 1):
        outfile4.write(str(i))
        for j in range(len(final_mu)):
            outfile4.write('\t' + str(round(final_mu[j][i-start], 5)))
        outfile4.write('\n')
    outfile4.close()

    # File 5 - responsibilities
    outfile_name5 = run_dir + 'Responsibilities.txt'
    outfile5 = open(outfile_name5, 'w')
    outfile5.write('Number\t')
    for k in range(K):
        k += 1
        outfile5.write('Cluster_' + str(k) + '\t')
    outfile5.write('N\tBit_vector\n')
    index_num = 1
    for bit_vect in X.n_occur:
        bv = ''.join(bit_vect)
        abundance = str(X.n_occur[bit_vect])
        outfile5.write(str(index_num) + '\t')
        for k in range(K):
            outfile5.write(str(round(resps[index_num-1][k], 3)) + '\t')
        outfile5.write(abundance + '\t' + bv + '\n')
        outfile5.write('\n')
        index_num += 1
    outfile5.close()

    # File 6 - Cluster proportions
    outfile_name6 = run_dir + 'Proportions.txt'
    outfile6 = open(outfile_name6, 'w')
    outfile6.write('Cluster, Obs Pi, Real pi \n')
    for k in range(K):
        obs_prob = str(round(final_obs_pi[k], 2))
        real_prob = str(round(final_real_pi[k], 2))
        outfile6.write(str(k+1) + ',' + obs_prob + ',' + real_prob + '\n')
    outfile6.close()

    # File 7 - Cluster proportions excluding all-zero bit vectors
    outfile_name7 = run_dir + 'Proportions_Excluding_Zeros.txt'
    outfile7 = open(outfile_name7, 'w')
    outfile7.write('Cluster, Proportion\n')
    
    total_resp_by_cluster = [0] * K
    total_weight = 0
    
    index_num = 1
    for bit_vect in X.n_occur:
        bv = ''.join(bit_vect)
        abundance = X.n_occur[bit_vect]
        
        if all(b == '0' for b in bit_vect): # Exclude all-zero bit vectors
            index_num += 1
            continue
            
        for k in range(K):
            total_resp_by_cluster[k] += resps[index_num-1][k] * abundance
        total_weight += abundance
        index_num += 1
    
    for k in range(K):
        proportion = total_resp_by_cluster[k] / total_weight if total_weight > 0 else 0
        outfile7.write(f"{k+1}, {round(proportion, 4)}\n")
    outfile7.close()

    # Plot 1 - log likelihood vs iteration number
    loglike_trace = go.Scatter(
        x=[(i+1) for i in range(len(log_like_list))],
        y=log_like_list,
        mode='lines'
    )
    loglike_layout = dict(xaxis=dict(title='Iteration'),
                          yaxis=dict(title='Log likelihood'))
    loglike_data = [loglike_trace]
    loglike_fig = dict(data=loglike_data, layout=loglike_layout)
    plotly.offline.plot(loglike_fig, filename=run_dir +
                        'LogLikes_Iterations.html',
                        auto_open=False)

    # Plot 2 - DMS mod rate for each base in each cluster
    DMSModRate_cluster_data = []
    xaxis_coords = [i for i in range(start, end+1)]
    for k in range(K):
        obs_prob = round(final_obs_pi[k], 2)
        real_prob = round(final_real_pi[k], 2)
        c_name = 'Cluster ' + str(k + 1) + ', obs p=' + str(obs_prob) + \
                 ', real p=' + str(real_prob)
        trace = go.Scatter(
            x=xaxis_coords,
            y=final_mu[k],
            name=c_name,
            mode='lines+markers'
        )
        DMSModRate_cluster_data.append(trace)
    DMSModRate_cluster_layout = dict(xaxis=dict(title='Position (BP)'),
                                     yaxis=dict(title='DMS mod rate'))
    DMSModRate_cluster_fig = dict(data=DMSModRate_cluster_data,
                                  layout=DMSModRate_cluster_layout)
    plotly.offline.plot(DMSModRate_cluster_fig, filename=run_dir +
                        'DMSModRate.html',
                        auto_open=False)

    # Plot 3 - Same as Plot 2, but in subplots
    cmap = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue'}  # Color map
    colors = [cmap[seq[i]] for i in range(len(seq))]
    ref_bases = [seq[i] for i in range(len(seq))]
    titles = ['Cluster ' + str(k+1) for k in range(K)]
    fig3 = tools.make_subplots(rows=K, cols=1, subplot_titles=titles)
    for k in range(K):
        trace = go.Bar(
            x=xaxis_coords,
            y=final_mu[k],
            text=ref_bases,
            marker=dict(color=colors),
            showlegend=False
        )
        fig3.append_trace(trace, k + 1, 1)
    plotly.offline.plot(fig3, filename=run_dir +
                        'DMSModRate_Clusters.html', auto_open=False)


def NumReads_File(sample_name, X, outplots_dir):
    """
    """
    outfile_name = outplots_dir + 'BitVectors_Filter.txt'
    outfile = open(outfile_name, 'w')
    outfile.write('Number of bit vectors used: ' + str(X.n_bitvectors) + '\n')
    outfile.write('Number of unique bit vectors used: ' +
                  str(X.n_unique_bitvectors) + '\n')
    outfile.write('Number of bit vectors discarded: ' +
                  str(X.n_discard) + '\n')
    outfile.close()


def LogLikes_File(sample_name, K, RUNS, log_likes, BICs,
                  best_run, outplots_dir):
    """
    """
    K_dir = outplots_dir + '/K_' + str(K) + '/'
    loglikes_file_name = K_dir + 'log_likelihoods.txt'
    loglikes_file = open(loglikes_file_name, 'w')
    loglikes_file.write('Run\tLog_likelihood\tBIC_score\n')
    for run in range(1, RUNS + 1):
        log_like = str(round(log_likes[run - 1], 2))
        BIC = str(round(BICs[run - 1], 2))
        line = str(run) + '\t' + log_like + '\t' + BIC + '\n'
        if run == best_run:
            line = str(run)+'-best' + '\t' + log_like + '\t' + BIC + '\n'
        loglikes_file.write(line)
    loglikes_file.close()


def Log_File(sample_name, NUM_RUNS, MIN_ITS,
             CONV_CUTOFF, INFO_THRESH, SIG_THRESH, exc_AC,
             norm_bases, K, time_taken, outplots_dir):
    """
    """
    now = datetime.datetime.now()
    log_file_name = outplots_dir + 'log.txt'
    log_file = open(log_file_name, 'w')
    log_file.write('Sample: ' + sample_name + '\n')
    log_file.write('Number of EM runs: ' + str(NUM_RUNS) + '\n')
    log_file.write('Minimum number of iterations: ' + str(MIN_ITS) + '\n')
    log_file.write('Convergence cutoff: ' + str(CONV_CUTOFF) + '\n')
    log_file.write('Informative bits threshold: ' + str(INFO_THRESH) + '\n')
    log_file.write('Signal threshold: ' + str(SIG_THRESH) + '\n')
    log_file.write('Exclude As and Cs?: ' + exc_AC + '\n')
    log_file.write('Num bases for normalization: ' + str(norm_bases) + '\n')
    log_file.write('Predicted number of clusters: ' + str(K) + '\n')
    log_file.write('Time taken: ' + str(time_taken) + ' mins\n')
    log_file.write('Finished at: ' + now.strftime("%Y-%m-%d %H:%M") + '\n')
    log_file.close()
