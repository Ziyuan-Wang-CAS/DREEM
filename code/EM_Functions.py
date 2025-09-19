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
from math import log
import numpy as np
import scipy.stats
import pandas as pd
import sys
sys.setrecursionlimit(10000)

def mu_der(mu_k, x_bar_k):
    """
    """
    epsilon = 1e-9
    mu_k_rev = mu_k[::-1]
    denom_k = calc_denom(0, mu_k, {}, {})
    denom_k_epsilon = denom_k[0] + epsilon
    denom_k_rev = calc_denom(0, mu_k_rev, {}, {})
    upd_mu = [(mu_k[i] * denom_k[1][i] * denom_k_rev[1][len(mu_k) - i - 1] /
              denom_k_epsilon) - x_bar_k[i] for i in range(len(mu_k))]
    return np.array(upd_mu)


def calc_denom(i, mu, denom_probs, s2_probs):
    """
        """
    if i in denom_probs:  # Already encountered
        return (denom_probs[i], s2_probs)
    elif i >= len(mu):  # Base case
        return (1, s2_probs)
    else:  # Make the calc
        s1 = calc_denom(i + 1, mu, denom_probs, s2_probs)[0]
        
        mu_next = np.clip(mu[i + 1: i + 4], 1e-6, 1 - 1e-6)
        s2 = (1.0 - mu_next).prod() * \
            calc_denom(i + 4, mu, denom_probs, s2_probs)[0]

        mu_i_clipped = np.clip(mu[i], 1e-6, 1 - 1e-6)
        denom_probs[i] = ((1 - mu_i_clipped) * s1) + (mu_i_clipped * s2)
        s2_probs[i] = s2
        return (denom_probs[i], s2_probs)


def is_distmuts_valid(bs):
    """
    """
    for i in range(len(bs)):
        if bs[i] == '1':
            try:
                if i - latest_mutbit_index < 4:
                    return False
            except NameError:  # This happens the first time we see a '1'
                None
            latest_mutbit_index = i
    return True


def is_surmuts_valid(bs):
    """
    """
    invalid_set = ['.1', '?1', '1.', '1?']
    for i in range(len(bs)):
        if bs[i:i + 2] in invalid_set:
            return False
    return True


def calc_nmuts_thresh(bv_filename):
    """
    """
    n_muts = pd.read_csv(bv_filename, sep='\t', skiprows=2,
                         usecols=['N_Mutations'], index_col=False)
    n_muts = n_muts['N_Mutations']
    mad = abs(n_muts - n_muts.median()).median()
    nmuts_thresh = n_muts.median() + (3 * mad / 0.6745)
    return int(round(nmuts_thresh))


def logpmf_function1(X, mu, ind, k):
    """
    """
    start, end = ind[0], ind[1]
    return scipy.stats.bernoulli.logpmf(X.BV_Matrix[start:end + 1], mu[k])


def logpmf_function2(log_pmf, denom, ind, k):
    """
    """
    start, end = ind[0], ind[1]
    return log_pmf[start:end + 1, k] - log(denom[k][0])


def calc_matrixIndices(N, K, cpus):
    """
    """
    calcsPerCPU = max(round(N * K / cpus), 1)
    inds, start = [], 0
    while start < N:
        coord = (start, start + calcsPerCPU - 1)
        inds.append(coord)
        start = start + calcsPerCPU
    return inds


def calc_BIC(N, D, K, log_like):
    """
    """
    return log(N) * D * K - (2 * log_like)
