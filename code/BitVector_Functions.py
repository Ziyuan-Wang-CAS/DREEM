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

Contains functions and a class used by Bit_Vector.py
"""


def Calc_Ambig_Reads(ref_seq, i, length, num_surBases):
    """
    Determines whether a deletion is ambiguous or not by looking at the
    sequence surrounding the deletion. Edge cases not handled right now.
    Args:
        ref_seq (string): Reference sequence
        i (int): 3' index of del at ref sequence
        length (int): Length of deletion
        num_surBases (int): Number of surrounding bases to consider
    Returns:
        boolean: Whether deletion is ambiguous or not
    """
    orig_del_start = i - length + 1
    orig_sur_start = orig_del_start - num_surBases
    orig_sur_end = i + num_surBases
    orig_sur_seq = ref_seq[orig_sur_start - 1: orig_del_start - 1] + \
        ref_seq[i:orig_sur_end]
    for new_del_end in range(i - length, i + length + 1):  # Alt del end points
        if new_del_end == i:  # Orig end point
            continue
        new_del_start = new_del_end - length + 1
        sur_seq = ref_seq[orig_sur_start - 1: new_del_start - 1] + \
            ref_seq[new_del_end:orig_sur_end]
        if sur_seq == orig_sur_seq:
            return True
    return False


def Parse_CIGAR(cigar_string):
    """
    Parse a CIGAR string
    Args:
        cigar_string (string): CIGAR string
    Returns:
        ops (list): List of operations. Each op is of type
        (length, description) such as ('37', 'M'), ('10', 'I'), (24, 'D'), etc.
    """
    import re
    ops = re.findall(r'(\d+)([A-Z]{1})', cigar_string)
    return ops


def Parse_PhredFile(qscore_filename):
    """
    Parse a file containing Phred Q Score info
    Args:
        qscore_filename (string): Path to Q Score file
    Returns:
        phred_qscore (dict): Mapping of ASCII symbol to Phred Q Score
    """
    phred_qscore = {}
    qscore_file = open(qscore_filename)
    qscore_file.readline()  # Ignore header line
    for line in qscore_file:
        line = line.strip().split()
        score, symbol = int(line[0]), line[1]
        phred_qscore[symbol] = score
    qscore_file.close()
    return phred_qscore


def Parse_FastaFile(fasta_file):
    """
    Parse a FASTA file
    Args:
        fasta_file (string): Path to FASTA file
    Returns:
        refs_seq (dict): Sequences of the ref genomes in the file
    """
    from Bio import SeqIO
    refs_seq = {}
    with open(fasta_file, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            refs_seq[record.id] = str(record.seq)
    return refs_seq


def Create_FastaFile(filename, ref_name, seq):
    """
    Create a file in FASTA format that's used by BLAST
    Args:
        filename (string): Output file name
        ref_name (string): Name of ref
        seq (string): Sequence to write to file
    """
    outfile = open(filename, 'w')
    outfile.write('>' + ref_name + '\n')
    outfile.write(seq + '\n')
    outfile.close()


class Mate():
    """
    Attributes of a mate in a read pair. Be careful about the order!
    """
    def __init__(self, line_split):
        self.QNAME = line_split[0]
        self.FLAG = line_split[1]
        self.RNAME = line_split[2]
        self.POS = int(line_split[3])
        self.MAPQ = int(line_split[4])
        self.CIGAR = line_split[5]
        self.RNEXT = line_split[6]
        self.PNEXT = int(line_split[7])
        self.TLEN = line_split[8]
        self.SEQ = line_split[9]
        self.QUAL = line_split[10]
        self.MDSTRING = line_split[11].split(":")[2]

    def __repr__(self):
        return self.QNAME + "-" + self.RNAME + "-" + str(self.POS)+"-" + \
            self.CIGAR+"-"+self.MDSTRING
