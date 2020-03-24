"""
Filter a fastq file
"""

import os
import sys
import argparse
from roblib import bcolors, stream_fastq
__author__ = 'Rob Edwards'



def filter_len(fq, length, verbose=False):
    """
    Filter by length.
    :param fq: array of seqs
    :param length: minimum length
    :param verbose: more output
    """
    if args.v:
        sys.stderr.write(f"{bcolors.GREEN}Filtering on length{bcolors.ENDC}\n")

    fqnew = []
    for s in fq:
        if len(s[2]) > length:
            fqnew.append(s)
    return fqnew


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='fastq file', required=True)
    parser.add_argument('-m', help='filter based on sequence length. Supply minimum length', type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    # just read the whole file into an array, and then we 
    #  can run serial filters
    fq = []
    for seqid, header, seq, scores in stream_fastq(args.f):
        fq.append([seqid, header, seq, scores])

    if args.m:
        fq = filter_len(fq, args.m, args.v)

    for s in fq:
        print(f"@{s[1]}\n{s[2]}\n+\n{s[3]}")


