"""
Count the kmers in a file and report entropy and eveness
"""

import os
import sys
import argparse

from roblib import bcolors
from itertools import product

from math import log2
from roblib import stream_fasta, rc, bcolors

def count_kmers(faf, k, verbose=False):
    """
    Count the kmers
    :param faf: fasta file
    :param k: kmer size
    :param verbose: more output
    :return: a dict of kmers
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Counting kmers (k={k}) in {faf}\n")

    kmers = {"".join(x) : 0 for x in product("ATGC", repeat=k)}
    for id, seq in stream_fasta(faf):
        rcseq = rc(seq)
        for x in kmers.keys():
            kmers[x] += seq.upper().count(x)
            kmers[x] += rcseq.upper().count(x)

    if verbose:
        sys.stderr.write(f"{bcolors.BLUE}\tDone counting kmers (k={k}) in {faf}\n")

    return kmers

def shannon(kmers, verbose=False):
    """
    Calculate the shannon entropy
    :param kmers: the kmer dictionary
    :param verbose: more output
    :return: the shannon entropy of the kmers
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Calculating Shannon's Entropy\n")
    t = sum(kmers.values())
    H = 0
    for x in kmers:
        H += (kmers[x] / t) * (log2(kmers[x]/t))

    return -H

def evenness(kmers, H=None, verbose=False):
    """
    Calculate the evenness
    :param kmers: the kmer dictionary
    :param H: shannon entropy (optional). If provided, we won't recalculate
    :param verbose: more output
    :return: the evenness of the kmers
    """

    if not H:
        H = shannon(kmers, verbose)
    S = len(kmers.keys())
    return H/log2(S)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count the kmers in a file and report entropy and eveness')
    parser.add_argument('-f', help='fasta file to count the entropy/evenness', required=True)
    parser.add_argument('-k', help='kmer size', required=True, type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    kmers = count_kmers(args.f, args.k, args.v)
    H = shannon(kmers, args.v)
    e = evenness(kmers, H, args.v)

    print(f"{args.f}\t{args.k}\t{H}\t{e}")