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
import numpy as np
import scipy.special
from scipy.optimize import newton_krylov
from multiprocessing.dummy import Pool as ThreadPool
import EM_Functions
import sys
sys.setrecursionlimit(10000)

def Run_EM(X, K, MIN_ITS, CONV_CUTOFF, CPUS):
    """
    Run the EM algorithm on the bit vector data contained in X.
    """
    conv_string = 'Log like converged after {:d} iterations'
    N, D = X.BV_Matrix.shape[0], X.BV_Matrix.shape[1]

    # Start and end coordinates of matrix for each thread
    calc_inds = EM_Functions.calc_matrixIndices(N, K, CPUS)
    
    avg_mut_rates = np.array([X.mut_popavg[d] for d in range(D)])
    mean_rate = np.mean(avg_mut_rates)
    perturbation_scale = mean_rate * 1.00  # Perturbation scale

    mu = np.asarray([avg_mut_rates for k in range(K)])
    
    # Initialize mutation rate vectors for each cluster
    for k in range(K):
        mu[k] = avg_mut_rates + np.random.normal(0, perturbation_scale, size=D)
        mu[k] = np.clip(mu[k], 1e-6, 1 - 1e-6)  # Ensure initialized mu is in valid range
    
    # Initialize cluster probabilities with a uniform distribution
    obs_pi = np.asarray([1.0 / K] * K)

    converged = False
    iteration = 1
    log_like_list, mu_list, obs_pi_list, real_pi_list = [], [], [], []

    momentum = 0.5  # Momentum for the EM algorithm
    prev_mu_update = np.zeros_like(mu)
    
    while not converged:  # Each iteration of the EM algorithm

        # Expectation step
        (resps, log_like, denom) = Exp_Step(X, K, mu, obs_pi, calc_inds, CPUS)

        # Maximization step with momentum
        (new_mu, obs_pi, real_pi) = Max_Step(X, K, mu, resps, denom)
        
        # Update mu with momentum
        mu_update = new_mu - mu
        mu = mu + mu_update + momentum * prev_mu_update
        prev_mu_update = mu_update
        mu = np.clip(mu, 0, 1) # Ensure mu is in valid range

        log_like_list.append(log_like)
        mu_list.append(mu)
        obs_pi_list.append(obs_pi)
        real_pi_list.append(real_pi)

        # Check if log like has converged
        if iteration >= MIN_ITS:  # At least min iterations has run
            prev_loglike = log_like_list[-2]
            diff = log_like - prev_loglike
            if diff <= CONV_CUTOFF:  # Converged
                converged = True
                print(conv_string.format(iteration))

        iteration += 1

    final_mu = mu_list[-1]
    final_obs_pi, final_real_pi = obs_pi_list[-1], real_pi_list[-1]

    BIC = EM_Functions.calc_BIC(N, D, K, log_like_list[-1])
    return (log_like_list, final_mu, final_obs_pi, final_real_pi, resps, BIC)


def Exp_Step(X, K, mu, pi, calc_inds, CPUS):
    """
    """
    N, D = X.BV_Matrix.shape[0], X.BV_Matrix.shape[1]
    log_pi = np.log(pi)
    log_pmf = np.zeros((N, D, K))
    denom = [EM_Functions.calc_denom(0, mu[k], {}, {}) for k in range(K)]

    input_array1 = [[X, mu, ind, k] for ind in calc_inds for k in range(K)]
    pool1 = ThreadPool(CPUS)
    logpmf_results1 = pool1.starmap(EM_Functions.logpmf_function1,
                                    input_array1)
    pool1.close()
    pool1.join()
    for i in range(len(logpmf_results1)):
        ind, k = input_array1[i][2], input_array1[i][3]
        start, end = ind[0], ind[1]
        log_pmf[start:end + 1, :, k] = logpmf_results1[i]

    min_log_pmf = -20  # Avoid underflow
    log_pmf = np.maximum(log_pmf, min_log_pmf)

    log_pmf = np.sum(log_pmf, axis=1)  # Sum of log - like taking product

    input_array2 = [[log_pmf, denom, ind, k] for ind in calc_inds
                    for k in range(K)]
    pool2 = ThreadPool(CPUS)
    logpmf_results2 = pool2.starmap(EM_Functions.logpmf_function2,
                                    input_array2)
    pool2.close()
    pool2.join()
    for i in range(len(logpmf_results2)):
        ind, k = input_array2[i][2], input_array2[i][3]
        start, end = ind[0], ind[1]
        log_pmf[start:end + 1, k] = logpmf_results2[i]

    log_resps_numer = np.add(log_pi, log_pmf)
    log_resps_denom = scipy.special.logsumexp(log_resps_numer, axis=1)
    log_resps = np.subtract(log_resps_numer.T, log_resps_denom).T
    resps = np.exp(log_resps)

    log_like = np.dot(log_resps_denom, X.BV_Abundance)
    return (resps, log_like, denom)


def Max_Step(X, K, mu, resps, denom):
    """
    """
    D = X.BV_Matrix.shape[1]
    mu, obs_pi, real_pi = np.zeros((K, D)), np.zeros(K), np.zeros(K)
    for k in range(K):
        N_k = np.sum(resps[:, k] * X.BV_Abundance)
        x_bar_k = np.sum((resps[:, k] * X.BV_Abundance *
                          X.BV_Matrix.T).T, axis=0) / N_k
        upd_mu = newton_krylov(lambda mu_k: EM_Functions.mu_der(mu_k, x_bar_k),
                               mu[k])
        mu[k] = upd_mu  # Mu with denom correction
        obs_pi[k] = N_k / X.n_bitvectors
    real_pi = [obs_pi[k] / denom[k][0] for k in range(K)]
    real_pi = real_pi / np.sum(real_pi)
    return (mu, obs_pi, real_pi)