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

Fold RNA sequence with RNAstructure using DMS contraints.
"""
import os
import pandas as pd
import numpy as np
from statistics import median as med
import BitVector_Functions
import json


def ConstraintFoldDraw(input_dir, clustMuFile, expUp, expDown, norm_bases):
    """
    """
    sample = clustMuFile.split('/')[-4]
    file_dir = os.path.dirname(clustMuFile)

    # Create trimmed reference fasta file
    clustMuFileContents = open(clustMuFile)
    first_line = clustMuFileContents.readline().strip()
    second_line = clustMuFileContents.readline().strip()
    clustMuFileContents.close()
    first_line_split = first_line.strip().split()
    ref_info = first_line_split[1]
    ref_file, ref = ref_info.split(';')[0], ref_info.split(';')[1]
    ref_file = input_dir + ref_file + '.fasta'
    refs_seq = BitVector_Functions.Parse_FastaFile(ref_file)
    entire_seq = refs_seq[ref]
    second_line_split = second_line.strip().split()
    indices = second_line_split[1].split(':')[0]
    start, end = int(indices.split(',')[0]), int(indices.split(',')[1])
    mod_start = max(1, start - expUp)
    mod_end = min(end + expDown, len(entire_seq))
    trim_seq = entire_seq[mod_start - 1:mod_end]
    trimref_filename = file_dir + '/' + ref + '_trimUp_' + str(expUp) + \
        '_trimDown_' + str(expDown) + '.fa'
    ref_name = ref + '_' + str(mod_start) + '_' + str(mod_end)
    BitVector_Functions.Create_FastaFile(trimref_filename, ref_name, trim_seq)

    # Gather mus for every k and normalize them
    clusts_mus = pd.read_csv(clustMuFile, sep='\t', skiprows=2,
                             index_col=False)
    rows, K = len(clusts_mus), len(clusts_mus.columns) - 1
    norm_clusts_mus = np.zeros((rows, K))
    for k in range(K):
        mus = clusts_mus['Cluster_' + str(k + 1)]
        norm_value = med(np.sort(mus)[-1:-(norm_bases+1):-1])  # Median of mus
        norm_mus = mus / norm_value  # Normalize the mus
        norm_mus[norm_mus > 1.0] = 1.0  # Cap at 1
        norm_clusts_mus[:, k] = norm_mus
    norm_clusts_mus = np.around(norm_clusts_mus, decimals=3)

    # Save normalized reactivity to json file for each cluster
    if 'Position' in clusts_mus.columns:
        positions = clusts_mus['Position'].tolist()
        for k in range(K):
            norm_data = []
            raw_mus = clusts_mus['Cluster_' + str(k + 1)].tolist()
            filtered_mus = np.copy(norm_clusts_mus[:, k])
            for i in range(rows):
                if isinstance(filtered_mus[i], str) and filtered_mus[i] == 'nan':
                    filtered_mus[i] = 0
                if isinstance(raw_mus[i], str) and raw_mus[i] == 'nan':
                    raw_value = 0
                else:
                    raw_value = raw_mus[i]
                entry = {
                    "position": positions[i],
                    "Raw_Reactivity": float(raw_value),
                    "Norm_Reactivity": float(filtered_mus[i]),
                }
                norm_data.append(entry)
            json_filename = file_dir + '/Reactivity_' + str(start) + '_' + str(end) + '_Cluster_' + str(k + 1) + '.json'
            with open(json_filename, 'w') as f:
                json.dump(norm_data, f, indent=4)

    # Drawing of structure for each k
    for k in range(K):
        clust_name = file_dir + '/' + sample + '-K' + str(K) + '_Cluster' + str(k+1)

        const_filename = clust_name + '_expUp_' + str(expUp) + '_expDown_' + \
            str(expDown) + '_const.txt'
        const_file = open(const_filename, 'w')
        for i in range(len(clusts_mus)):
            pos = clusts_mus['Position'][i]
            mod_pos = pos - mod_start + 1  # Pos wrt trimmed seq
            mu = str(norm_clusts_mus[i][k])
            if mu == 'nan':  # Happens in UT
                mu = '0'
            if entire_seq[pos-1] == 'A' or entire_seq[pos-1] == 'C':
                mu = '-999'
            if mod_pos > 0 and mod_start <= pos <= mod_end:  # Can be < 0 because its wrt trimmed seq
                const_file.write(str(mod_pos) + '\t' + mu + '\n')
        const_file.close()

        # Folding using RNAstructure
        ct_filename = clust_name + '_expUp_' + str(expUp) + '_expDown_' + \
            str(expDown) + '.ct'
        dot_filename = clust_name + '_expUp_' + str(expUp) + '_expDown_' + \
            str(expDown) + '.dot'
        pic_filename = clust_name + '_basesExpanded_' + str(expUp) + '.ps'
        fold_command = 'Fold -m 3 ' + trimref_filename + ' -dms ' + \
            const_filename + ' ' + ct_filename + ' -md 350'
        ct2dot_command = 'ct2dot ' + ct_filename + ' ALL ' + \
            dot_filename
        draw_command = 'draw ' + dot_filename + ' ' + \
            pic_filename + ' -S ' + const_filename
        os.system(fold_command)
        os.system(ct2dot_command)
        os.system(draw_command)

        # Delete unnecessary files
        os.system('rm ' + const_filename)
        # os.system('rm ' + ct_filename)
        # os.system('rm ' + dot_filename)

    os.system('rm ' + trimref_filename)
