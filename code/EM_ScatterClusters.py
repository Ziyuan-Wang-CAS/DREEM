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
import os
import pandas as pd
import numpy as np
import plotly
import plotly.graph_objs as go
from scipy.stats import linregress
import BitVector_Functions


def Scatter_Clusters(input_dir, clustMuFile):
    """
    """
    file_dir = os.path.dirname(clustMuFile)
    clustMuFileContents = open(clustMuFile)
    first_line = clustMuFileContents.readline().strip()
    clustMuFileContents.close()
    first_line_split = first_line.strip().split()
    ref_info = first_line_split[1]
    ref_file, ref = ref_info.split(';')[0], ref_info.split(';')[1]
    ref_file = input_dir + ref_file + '.fasta'
    refs_seq = BitVector_Functions.Parse_FastaFile(ref_file)
    entire_seq = refs_seq[ref]
    clusts_mus = pd.read_csv(clustMuFile, sep='\t', skiprows=2,
                             index_col=False)
    positions = clusts_mus['Position']
    valid_indices1 = [i for i in range(len(positions)) if
                      entire_seq[positions[i] - 1] == 'A' or
                      entire_seq[positions[i] - 1] == 'C' or
                      entire_seq[positions[i] - 1] == 'T' or
                      entire_seq[positions[i] - 1] == 'G'] # Include Ts and Gs
    K = len(clusts_mus.columns) - 1
    r2_filename = file_dir + '/NormMu_RSquare.txt'
    r2_file = open(r2_filename, 'w')
    r2_file.write('ClusterA-ClusterB:R-Square,P-value\n')
    for k1 in range(K-1):
        for k2 in range(k1 + 1, K):
            mus1 = clusts_mus['Cluster_' + str(k1 + 1)]
            mus2 = clusts_mus['Cluster_' + str(k2 + 1)]
            i1 = np.where(mus1 > 0)[0]  # Non zero indices
            i2 = np.where(mus2 > 0)[0]  # Non zero indices
            valid_indices2 = [i for i in i1 if i in i2]  # Present in both

            mus1 = [mus1[i] for i in range(len(mus1)) if i in valid_indices1 and
                    i in valid_indices2]
            mus2 = [mus2[i] for i in range(len(mus2)) if i in valid_indices1 and
                    i in valid_indices2]
            if len(mus1) == 0 or len(mus2) == 0:  # Like in UT
                continue

            norm_value1 = np.sort(mus1)[-1]
            norm_mus1 = mus1 / norm_value1
            norm_value2 = np.sort(mus2)[-1]
            norm_mus2 = mus2 / norm_value2

            # Plot 1 - normalized mus
            slope1, intercept1, r1, p1, std_err1 = linregress(norm_mus1,
                                                              norm_mus2)
            r1 = r1**2
            r1, p1 = round(r1, 2), round(p1, 2)
            trace1 = go.Scatter(
                x=norm_mus1,
                y=norm_mus2,
                mode='markers',
                marker=dict(color='rgb(49,130,189)'),
                showlegend=False
            )
            trace2 = go.Scatter(
                x=[0, 1],
                y=[0, 1],
                mode='lines',
                line=dict(color='rgb(49,130,189)'),
                showlegend=False
            )
            title1 = 'R-squared: ' + str(r1) + ', p-value: ' + str(p1)
            layout1 = dict(title=title1,
                           xaxis=dict(title='Cluster '+str(k1 + 1)),
                           yaxis=dict(title='Cluster '+str(k2 + 1)))
            data1 = [trace1, trace2]
            fig1 = dict(data=data1, layout=layout1)
            fname1 = file_dir + '/Cluster' + str(k1 + 1) + '_Cluster' + \
                str(k2 + 1) + '_normmus.html'
            plotly.offline.plot(fig1, filename=fname1, auto_open=False)

            # Plot 2 - mus
            slope2, intercept2, r2, p2, std_err2 = linregress(mus1, mus2)
            r2 = r2 ** 2
            r2, p2 = round(r2, 2), round(p2, 2)
            m1 = max(max(mus1), max(mus2))
            trace1_1 = go.Scatter(
                x=mus1,
                y=mus2,
                mode='markers',
                marker=dict(color='rgb(49,130,189)'),
                showlegend=False
            )
            trace2_1 = go.Scatter(
                x=[0, m1],
                y=[0, m1],
                mode='lines',
                line=dict(color='rgb(49,130,189)'),
                showlegend=False
            )
            title2 = 'R-squared: ' + str(r2) + ', p-value: ' + str(p2)
            layout2 = dict(title=title2,
                           xaxis=dict(title='Cluster '+str(k1 + 1)),
                           yaxis=dict(title='Cluster '+str(k2 + 1)))
            data2 = [trace1_1, trace2_1]
            fig2 = dict(data=data2, layout=layout2)
            fname2 = file_dir + '/Cluster' + str(k1 + 1) + '_Cluster' + \
                str(k2 + 1) + '_mus.html'
            plotly.offline.plot(fig2, filename=fname2, auto_open=False)

            r2_file.write('Cluster' + str(k1+1) + '-Cluster' + str(k2+1) +
                          ':' + str(r1) + ',' + str(p1) + '\n')
    r2_file.close()
