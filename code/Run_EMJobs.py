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
import EM_Files
from EM_Algorithm import Run_EM
import EM_Plots
import sys
sys.setrecursionlimit(10000)


def Run_EMJob(X, bvfile_basename, MIN_ITS, INFO_THRESH, CONV_CUTOFF,
              SIG_THRESH, outplot_dir, K, CPUS, run):

    if K == 1:
        EM_Plots.NumReads_File(bvfile_basename, X, outplot_dir)

    EM_res = Run_EM(X, K, MIN_ITS, CONV_CUTOFF, CPUS)
    log_like_list, final_mu, final_obs_pi, final_real_pi, resps, BIC = EM_res

    EM_Plots.Run_Plots(bvfile_basename, X, K, log_like_list, final_mu,
                       final_obs_pi, final_real_pi, resps, BIC, outplot_dir,
                       run)
