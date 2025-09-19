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

Class for a Bit Vector matrix
"""
import numpy as np
import random


class BV_Object():
    """
    """
    def __init__(self, bit_vectors, mut_popavg, n_discard, ref_file,
                 ref, seq, infiles_dir, indices):
        BV_Matrix, BV_Abundance, n_occur = [], [], {}
        for bit_vector in bit_vectors:
            bit_vector = tuple(bit_vector)  # Change to a tuple
            if bit_vector in n_occur:
                n_occur[bit_vector] += 1
            else:
                n_occur[bit_vector] = 1
        for bit_vector in n_occur:
            bv = np.array(list(map(float, bit_vector)))  # Convert to float
            BV_Matrix.append(bv)
            BV_Abundance.append(n_occur[bit_vector])

        BV_Matrix = np.array(BV_Matrix)
        BV_Abundance = np.array(BV_Abundance)
        self.BV_Matrix = BV_Matrix  # Only unique bit vectors
        self.BV_Abundance = BV_Abundance  # Abundance of each bit vector
        self.n_occur = n_occur
        self.n_bitvectors = len(bit_vectors)
        self.n_unique_bitvectors = len(n_occur.keys())
        self.n_discard = n_discard
        self.mut_popavg = mut_popavg
        self.ref = ref
        self.ref_file = ref_file
        self.seq = seq
        self.infiles_dir = infiles_dir
        self.indices = indices
