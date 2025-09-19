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

Create outputs of the bitvector step: text file, plots, etc.
"""
import datetime
import pandas as pd
import plotly
import plotly.graph_objs as go
from plotly import tools


def writeOutputFiles(sample_name, ref_name, num_reads, outfile_dir,
                     outplots_dir, refs_seq, start, end, mod_bases,
                     mut_bases, delmut_bases, info_bases, cov_bases,
                     ref_filename, sur_bases, qscore_filename,
                     qscore_cutoff, time_taken):
    """
    Write Bitvector info to output file
    Create all relevant plots
    Args:
        sample_name (string): Name of sample
        ref_name (string): Name of the reference genome
        num_reads (dict): Number of reads per ref
        outfile_dir (string): Path to output file directory
        outplots_dir (string): Path to output plots directory
        refs_seq (dict): Ref genome sequences
        start (int): Start position
        end (int): End position
        mod_bases (dict): Number of times a base was modified to at a pos
        mut_bases (dict): Number of times mut occurred at a pos
        delmut_bases (dict): Number of times mut and del occurred at a pos
        info_bases (dict): Number of data points with info at a pos
        cov_bases (dict): Number of times a base was covered
        ref_filename (string): Path to ref file
        sur_bases (int): Number of surrounding bases for deletion
        qscore_filename (string): Q score-symbol mapping file
        qscore_cutoff (int): Cutoff for valid base
        time_taken (float): Time taken for creating the bit vectors
    """
    now = datetime.datetime.now()
    bases = ['A', 'T', 'G', 'C', 'N']
    cmap = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue', 'N': 'gray'}  # Color map

    # Write to log file
    log_file_name = outplots_dir + sample_name + '_' + ref_name + '_' + \
        'log.txt'
    log_file = open(log_file_name, 'w')
    log_file.write('Sample: ' + sample_name + '\n')
    log_file.write('Reference file: ' + ref_filename + '\n')
    log_file.write('Reference genomes: ' + str(list(refs_seq.keys())) + '\n')
    log_file.write('Num surrounding bases for del: ' + str(sur_bases) + '\n')
    log_file.write('Q score file: ' + qscore_filename + '\n')
    log_file.write('Q score cutoff: ' + str(qscore_cutoff) + '\n')
    log_file.write('Time taken: ' + str(time_taken) + ' mins\n')
    log_file.write('Finished at: ' + now.strftime("%Y-%m-%d %H:%M") + '\n')
    log_file.close()

    # Write bitvector txt file and create plots for each ref
    for ref in refs_seq:  # Each ref genome
        file_base_name = sample_name + '_' + ref + '_' + str(start) + \
            '_' + str(end) + '_'
        ref_seq = refs_seq[ref]
        bv_filename = outfile_dir + file_base_name + 'bitvectors.txt'
        n_muts = pd.read_csv(bv_filename, sep='\t', skiprows=2,
                             usecols=['N_Mutations'], index_col=False)
        n_muts = n_muts['N_Mutations']

        # Plot 1 - Read coverage
        xaxis_coordinates = [i for i in range(start, end + 1)]
        read_cov = []
        for pos in range(start, end + 1):
            try:
                cov_frac = cov_bases[ref][pos] / num_reads[ref]
            except ZeroDivisionError:
                cov_frac = 0.0
            read_cov.append(cov_frac)
        cov_trace = go.Bar(
                x=xaxis_coordinates,
                y=read_cov
        )
        cov_data = [cov_trace]
        cov_layout = go.Layout(
            title='Read coverage: ' + sample_name +
            ', Number of bit vectors: ' + str(num_reads[ref]),
            xaxis=dict(title='Position'),
            yaxis=dict(title='Coverage fraction')
        )
        cov_fig = go.Figure(data=cov_data, layout=cov_layout)
        plotly.offline.plot(cov_fig, filename=outplots_dir + file_base_name +
                            'read_coverage.html', auto_open=False)

        # Plot 2 - Bar chart of abundance of modified bases
        modbases_data = []
        for base in bases:
            y_list = [mod_bases[ref][base][pos] for pos in
                      range(start, end + 1)]
            trace = go.Bar(
                    x=xaxis_coordinates,
                    y=y_list,
                    name=base,
                    marker=dict(color=cmap[base])
            )
            modbases_data.append(trace)
        modbases_layout = go.Layout(
            title='DMS modifications: ' + sample_name,
            xaxis=dict(title='Position'),
            yaxis=dict(title='Abundance'),
            barmode='stack'
        )
        modbases_fig = go.Figure(data=modbases_data, layout=modbases_layout)
        plotly.offline.plot(modbases_fig, filename=outplots_dir +
                            file_base_name + 'DMS_mutations.html',
                            auto_open=False)

        # Plot 3 - Histogram of number of mutations per read
        mut_hist_data = [go.Histogram(x=n_muts)]
        mut_hist_layout = go.Layout(
            title='Mutations: ' + sample_name,
            xaxis=dict(title='Number of mutations per read'),
            yaxis=dict(title='Abundance')
        )
        mut_hist_fig = go.Figure(data=mut_hist_data, layout=mut_hist_layout)
        plotly.offline.plot(mut_hist_fig, filename=outplots_dir +
                            file_base_name + 'mutation_histogram.html',
                            auto_open=False)

        # Plot 4 - Pop avg of mutated bases
        delmut_y, mut_y = [], []
        for pos in range(start, end + 1):
            try:
                delmut_frac = delmut_bases[ref][pos] / info_bases[ref][pos]
                mut_frac = mut_bases[ref][pos] / info_bases[ref][pos]
            except ZeroDivisionError:
                delmut_frac = 0.0
                mut_frac = 0.0
            delmut_y.append(delmut_frac)
            mut_y.append(mut_frac)
        colors = [cmap[ref_seq[i - 1]] for i in range(start, end + 1) if
                  i < len(ref_seq)]
        ref_bases = [ref_seq[i - 1] for i in range(start, end + 1) if
                  i < len(ref_seq)]
        delmut_trace = go.Bar(
                x=xaxis_coordinates,
                y=delmut_y,
                text=ref_bases,
                marker=dict(color=colors),
                showlegend=False
        )
        mut_trace = go.Bar(
                x=xaxis_coordinates,
                y=mut_y,
                text=ref_bases,
                marker=dict(color=colors),
                showlegend=False
        )
        title1 = 'Mismatches + Deletions: ' + sample_name
        title2 = 'Mismatches: ' + sample_name
        mut_fig = tools.make_subplots(rows=2, cols=1,
                                      subplot_titles=(title1, title2))
        mut_fig.append_trace(delmut_trace, 1, 1)
        mut_fig.append_trace(mut_trace, 2, 1)
        mut_fig['layout']['xaxis1'].update(title='Position')
        mut_fig['layout']['xaxis2'].update(title='Position')
        mut_fig['layout']['yaxis1'].update(title='Mutational fraction', range=[0, 0.1])
        mut_fig['layout']['yaxis2'].update(title='Mutational fraction', range=[0, 0.1])
        plotly.offline.plot(mut_fig, filename=outplots_dir + file_base_name +
                            'pop_avg.html', auto_open=False)

        '''
        # SIG_THRESH calculation for base 'T'
        T_mut = []
        for pos in range(start, end + 1):
            if pos in delmut_bases[ref] and ref_seq[pos - 1] == 'T':
                try:
                    T_mut.append(delmut_bases[ref][pos]/info_bases[ref][pos])
                except ZeroDivisionError:
                    T_mut.append(0)
        if T_mut:
            T_mut.sort(reverse=True) 
            index_90 = int(len(T_mut) * 0.9) 
            SIG_THRESH = round(T_mut[index_90], 4)
            print("Calculating signals threshold...done.\nSignals threshold: ", SIG_THRESH)
        '''
