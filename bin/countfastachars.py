#!/usr/bin/env python

"""
Count the characters in a fasta file. We summarize the longest and shortest reads
and the N50 of the data set.
"""

import os
import sys
import argparse
from roblib import read_fasta

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='fasta file', required=True)
    args = parser.parse_args()

    fa=read_fasta(args.f)
    lens=[len(fa[i]) for i in fa]
    lens.sort()
    length=sum(lens)
    n50 = None
    n75 = None
    n50name = ""
    n75name = ""
    shortestname = ""
    longestname = ""
    c=0

    for i in fa:
        if len(fa[i]) == lens[0]:
            shortestname = i
        if len(fa[i]) == lens[-1]:
            longestname = i
        c += len(fa[i])
        if not n50 and c>= length*0.5:
            n50 = len(fa[i])
            n50name = i
        if not n75 and c>=length*0.75:
            n75 = len(fa[i])
            n75name = i
    
    print("Total length: {}\nShortest: {} ({})\nLongest: {} ({})\n\nN50: {} ({})\nN75: {} ({})".format(
        length, lens[0], shortestname, lens[-1], longestname, n50, n50name, n75, n75name
    ))
