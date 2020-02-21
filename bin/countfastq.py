#!/usr/bin/env python3

"""
Count the characters in a fasta file. We summarize the longest and shortest reads
and the N50 of the data set.
"""

import os
import sys
import argparse
from roblib import stream_fastq

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-t', help='tab separated summmary of name, total len, shortest, longest, n50, n75', action="store_true")
    parser.add_argument('-s', help='(deprecated). Same as -t', action="store_true")
    args = parser.parse_args()

    lens = []
    for (sid, label, seq, qual) in stream_fastq(args.f):
        lens.append(len(seq))
    lens.sort()
    length=sum(lens)

    len_so_far = 0
    n50 = None
    n75 = None

    for i in lens:
        len_so_far += i
        if not n50 and len_so_far >= length * 0.5:
            n50 = i
        if not n75 and len_so_far >= length * 0.75:
            n75 = i

    if args.s or args.t:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(args.f, len(lens), length, lens[0], lens[-1], n50, n75))
    else:
        print("Number of sequences: {}\nTotal length: {}\nShortest: {}\nLongest: {}\nN50: {}\nN75: {}".format(
            len(lens), length, lens[0], lens[-1], n50, n75,
        ))
