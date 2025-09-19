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

Read the bitvector text file and create object X
Read the FASTA file to get the ref genome sequence
Filtering of the bit vectors is done here using various criteria
Changing of . and ? to 0 is done here
"""
import EM_Class
import EM_Functions
import numpy as np
import json
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import os
import time
import plotly
import plotly.graph_objs as go

fold_change_thresh = 1.5
q_value_thresh = 0.00001
epsilon = 1e-8 # Additive constant to avoid division by zero

def Load_BitVectors(bv_file, INFO_THRESH, SIG_THRESH, exc_AC, output_dir, ctrl):
    """
    """
    bases = ['A', 'T', 'G', 'C']
    bit_strings, mut_popavg, n_discard = [], {}, 0
    f, f1, f2, f3, f4 = 0, 0, 0, 0, 0

    bv_fileobj = open(bv_file)
    bvfile_contents = bv_fileobj.readlines()
    bv_fileobj.close()

    first_line = bvfile_contents[0]
    first_line_split = first_line.strip().split()
    ref_info, seq = first_line_split[1], first_line_split[2]
    ref_file, ref = ref_info.split(';')[0], ref_info.split(';')[1]

    second_line = bvfile_contents[1]
    second_line_split = second_line.strip().split()
    indices = second_line_split[1].split(':')[0]

    l = len(bvfile_contents[3].strip().split()[1])  # Len of 1st bit string
    nmuts_min = int(round(0.1 * l))
    nmuts_thresh = max(nmuts_min, EM_Functions.calc_nmuts_thresh(bv_file))
    print('Mutations threshold:', nmuts_thresh)

    for i in range(3, len(bvfile_contents)):
        f += 1
        line = bvfile_contents[i].strip().split()
        bit_string = line[1]
        n_mut = float(line[2])

        # Replace bases with 1
        for base in bases:
            bit_string = bit_string.replace(base, '1')

        # Filter 1 - Number of mutations
        if n_mut > nmuts_thresh:
            n_discard += 1
            f1 += 1
            continue

        # Filter 2 - Fraction of informative bits
        if (bit_string.count('.') + bit_string.count('?') +
           bit_string.count('N')) >= INFO_THRESH * len(bit_string):
            n_discard += 1
            f2 += 1
            continue

        # Filter 3 - Distance between mutations
        if not EM_Functions.is_distmuts_valid(bit_string):
            n_discard += 1
            f3 += 1
            continue

        # Filter 4 - Bits surrounding mutations
        if not EM_Functions.is_surmuts_valid(bit_string):
            n_discard += 1
            f4 += 1
            continue

        bit_strings.append(bit_string)

    print('Total bit vectors:', f)
    print('Bit vectors removed because of too many mutations: ', f1)
    print('Bit vectors removed because of too few informative bits: ', f2)
    print('Bit vectors removed because of mutations close by: ', f3)
    print('Bit vectors removed because of no info around mutations: ', f4)

    D = len(bit_strings[0])
    thresh_pos = []  # Positions above or below signal threshold
    for d in range(D):  # Each position of interest in the genome
        bits_list = [bs[d] for bs in bit_strings]  # List of bits at that pos
        noinfo_count = bits_list.count('.') + bits_list.count('?') + \
            bits_list.count('N')
        info_count = len(bits_list) - noinfo_count  # Num of informative bits
        try:
            mut_prob = bits_list.count('1') / info_count
        except ZeroDivisionError:
            mut_prob = 0
        mut_popavg[d] = mut_prob

    if ctrl == 'True':
        # ------------------- Recording of control group stats -------------------
        control_stats = {}
        for i in range(D):
            bits_list = [bs[i] for bs in bit_strings]
            noinfo_count = bits_list.count('.') + bits_list.count('?') + bits_list.count('N')
            valid_count = len(bits_list) - noinfo_count
            mut_count = bits_list.count('1')
            control_stats[i] = {"mut_count": mut_count, "total_count": valid_count}
        
        control_stats_file = output_dir + 'control_group_stats.json'
        with open(control_stats_file, 'w') as f:
            json.dump(control_stats, f, indent=4)
        print('Control group stats file saved to: ', control_stats_file)

        base_color = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue', 'N': 'gray'}
        start_pos = int(indices.split(',')[0])
        xaxis = [start_pos + i for i in range(D)]
        yaxis = [mut_popavg[pos] for pos in range(D)]
        ref_bases = list(seq.strip())[:D]
        colors = [base_color[base] for base in ref_bases]

        trace = go.Bar(
            x=xaxis,
            y=yaxis,
            text=ref_bases,
            marker=dict(color=colors),
            showlegend=False
        )
        layout = go.Layout(
            title="Control Group Mutation Population Average",
            xaxis=dict(title="Position"),
            yaxis=dict(title="Mutation Rate")
        )
        fig = go.Figure(data=[trace], layout=layout)
        plot_filename = output_dir + "control_mut_pop_avg.html"
        plotly.offline.plot(fig, filename=plot_filename, auto_open=False)
        print('Control group mut_popavg plot saved to: ', plot_filename)

    else:
        # ------------------- Screening of signals in experimental group -------------------
        control_stats_file = output_dir + 'control_group_stats.json'
        while not os.path.exists(control_stats_file):
            print('Warning: Control group stats file not found. Waiting for file import.')
            time.sleep(10)
        with open(control_stats_file, 'r') as f:
            control_group_stats = json.load(f)

        exp_stats = {}
        for i in range(D):
            bits_list = [bs[i] for bs in bit_strings]
            noinfo_count = bits_list.count('.') + bits_list.count('?') + bits_list.count('N')
            valid_count = len(bits_list) - noinfo_count
            mut_count = bits_list.count('1')
            exp_stats[i] = {"mut_count": mut_count, "total_count": valid_count}
        
        p_values = []
        positions_tested = list(range(D))
        for d in positions_tested:
            control_mut = control_group_stats[str(d)]["mut_count"]
            control_total = control_group_stats[str(d)]["total_count"]
            exp_mut = exp_stats[d]["mut_count"]
            exp_total = exp_stats[d]["total_count"]
            table = [[exp_mut, exp_total - exp_mut],
                     [control_mut, control_total - control_mut]]
            _, p = fisher_exact(table, alternative='greater')
            p_values.append(p)
        
        _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
        
        thresh_pos = []

        for d, q in enumerate(q_values):
            exp_total = exp_stats[d]["total_count"]
            ctrl_total = control_group_stats[str(d)]["total_count"]

            if exp_total == 0 or ctrl_total == 0:
                fold_change = 0
            else:
                exp_rate = exp_stats[d]["mut_count"] / exp_total
                ctrl_rate = control_group_stats[str(d)]["mut_count"] / ctrl_total

                fold_change = (exp_rate + epsilon) / (ctrl_rate + epsilon)
            
            if q >= q_value_thresh and fold_change < fold_change_thresh:
                mut_popavg[d] = 0
                thresh_pos.append(d)
        
        base_color = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue', 'N': 'gray'}
        start_pos = int(indices.split(',')[0])
        xaxis = [start_pos + i for i in range(D)]
        yaxis = [mut_popavg[pos] for pos in range(D)]

        ref_bases = list(seq.strip())[:D]
        colors = [base_color[base] for base in ref_bases]

        trace = go.Bar(
            x=xaxis,
            y=yaxis,
            text=ref_bases,
            marker=dict(color=colors),
            showlegend=False
        )
        layout = go.Layout(
            title="Experimental Group Mutation Population Average",
            xaxis=dict(title="Position"),
            yaxis=dict(title="Mutation Rate")
        )
        fig = go.Figure(data=[trace], layout=layout)
        plot_filename = output_dir + "experimental_mut_pop_avg.html"
        plotly.offline.plot(fig, filename=plot_filename, auto_open=False)
        print('Experimental group mut_popavg plot saved to: ', plot_filename)

    for i in range(len(bit_strings)):  # Change . and ? to 0, noise to 0
        bit_string = bit_strings[i]
        bit_string = bit_string.replace('?', '0')
        bit_string = bit_string.replace('.', '0')
        bit_string = bit_string.replace('N', '0')

        # Suppressing data from As and Cs
        if exc_AC == 'True':
            bit_string = list(bit_string)
            j = 0
            while j < len(bit_string):
                if seq[j] == 'A' or seq[j] == 'C':
                    bit_string[j] = '0'
                j += 1
            bit_string = ''.join(bit_string)

        bit_string = np.array(list(bit_string))
        bit_string[thresh_pos] = '0'
        bit_string = ''.join(bit_string)

        bit_strings[i] = bit_string

    # Update mutation population average for EM clustering
    mut_popavg = {}
    for d in range(D):
        bits_list = [bs[d] for bs in bit_strings]
        info_count = len(bits_list)
        try:
            mut_prob = bits_list.count('1') / info_count
        except ZeroDivisionError:
            mut_prob = 0
        mut_popavg[d] = mut_prob

    X = EM_Class.BV_Object(bit_strings, mut_popavg, n_discard, ref_file,
                           ref, seq, output_dir, indices)
    return X