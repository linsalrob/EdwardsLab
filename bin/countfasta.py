#!/usr/bin/env python3

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
    parser.add_argument('-f', nargs='+', help='fasta file', required=True)
    parser.add_argument('-l', help='list the lengths for each sequence (default = not to)', action='store_true')
    parser.add_argument('-m', help='minimum length fo be inclued', type=int, default=0)
    parser.add_argument('-t', help='tab separated output. Fields: [# seqs, total bp, shortest, longest, N50, N75]', action='store_true')
    args = parser.parse_args()


    for faf in args.f:
        fa=read_fasta(faf)

        if len(fa.keys()) == 1 and list(fa.keys())[0] == '':
            sys.stderr.write(f"No sequences found in {faf}\n")
            sys.exit(0)

        if args.l:
            for i in fa:
                print("{}\t{}".format(i, len(fa[i])))
            print()

        lensall=[len(fa[i]) for i in fa]
        lens = list(filter(lambda x: x > args.m, lensall))
        lens.sort()
        length=sum(lens)

        len_so_far = 0
        n50 = None
        n75 = None
        auN = 0
        for i in lens:
            len_so_far += i
            if not n50 and len_so_far >= length * 0.5:
                n50 = i
            if not n75 and len_so_far >= length * 0.75:
                n75 = i
            auN += i**2

        auN /= length

        if args.t:
            print("\t".join(map(str, [faf, len(lens), length, lens[0], lens[-1], n50, n75, auN])))
        else:
            print("File name: {}\nNumber of sequences: {}\nTotal length: {}\nShortest: {}\nLongest: {}\nN50: {}\nN75: {}\nauN: {}".format(
                faf, len(lens), length, lens[0], lens[-1], n50, n75, auN
            ))

