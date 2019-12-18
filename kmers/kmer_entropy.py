"""
Count the kmers in a file and report entropy and eveness
"""

import os
import sys
import argparse

from roblib import bcolors
from itertools import product

import json
from math import log2
from roblib import stream_fasta, rc, bcolors, stream_fastq

def count_kmers(faf, type, k, jsonout=None, verbose=False):
    """
    Count the kmers
    :param faf: fasta file
    :param type: str either fasta or fastq
    :param k: kmer size
    :param verbose: more output
    :return: a dict of kmers
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Counting kmers (k={k}) in {faf}\n")

    kmers = {}

    if type == "fasta":
        for id, seq in stream_fasta(faf):
            rcseq = rc(seq)
            posn = 0
            while posn < len(seq) - k - 1:
                kmers[seq[posn:posn+k]] = kmers.get(seq[posn:posn+k], 0) + 1
                kmers[rcseq[posn:posn + k]] = kmers.get(rcseq[posn:posn + k], 0) + 1
                posn += 1

    if type == "fastq":
        for id, fullid, seq, qual in stream_fastq(faf):
            rcseq = rc(seq)
            posn = 0
            while posn < len(seq) - k - 1:
                kmers[seq[posn:posn+k]] = kmers.get(seq[posn:posn+k], 0) + 1
                kmers[rcseq[posn:posn + k]] = kmers.get(rcseq[posn:posn + k], 0) + 1
                posn += 1

    if jsonout:
        if verbose:
            sys.stderr.write(f"{bcolors.BLUE}\tWriting to {jsonout}\n")
        with open(jsonout, 'w') as out:
            json.dump({faf : kmers}, out)

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
    parser.add_argument('-f', help='fasta file to count the entropy/evenness')
    parser.add_argument('-q', help='fastq file to count the entropy/evenness')
    parser.add_argument('-k', help='kmer size', required=True, type=int)
    parser.add_argument('-j', help='json output for kmer counts')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.f:
        kmers = count_kmers(args.f, 'fasta', args.k, args.j, args.v)
    elif args.q:
        kmers = count_kmers(args.q, 'fastq', args.k, args.j, args.v)
    else:
        sys.stderr.write(f"{bcolors.RED}FATAL: Please supply either a fasta file or a fastq file{bcolors.ENDC}\n")
        sys.exit(-1)

    H = shannon(kmers, args.v)
    e = evenness(kmers, H, args.v)

    print(f"{args.f}\t{args.k}\t{H}\t{e}")