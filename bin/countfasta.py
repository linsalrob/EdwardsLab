#!/usr/bin/env python3

"""
Count the characters in a fasta file. We summarize the longest and shortest reads
and the N50 of the data set.
"""

import os
import sys
import argparse
from roblib import read_fasta, bcolors

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', nargs='+', help='fasta file')
    parser.add_argument('-d', nargs='+', help='directory of fasta files')
    parser.add_argument('-l', help='list the lengths for each sequence (default = not to)', action='store_true')
    parser.add_argument('-m', help='minimum length fo be inclued', type=int, default=0)
    parser.add_argument('-n', help='do NOT print the summary at the end', action='store_true')
    parser.add_argument('-t', help='tab separated output. Fields: [SeqID, # seqs, total bp, shortest, longest, N50, N75, auN]', action='store_true')
    parser.add_argument('-v', help='more output', action='store_true')
    args = parser.parse_args()

    if not args.f and not args.d:
        sys.stderr.write(f"{bcolors.RED}FATAL: Please specify either -f or -d or use -h for more help{bcolors.ENDC}\n")
        sys.exit(1)

    if args.f:
        files = args.f
    else:
        files = []

    if args.d:
        for subdir in args.d:
            for f in os.listdir(subdir):
                if not (f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fna')):
                    sys.stderr.write(f"Skipped {f}: does not look like a fasta file\n")
                    continue
                files.append(os.path.join(subdir, f))
                if args.v:
                    sys.stderr.write(f"Added {files[-1]}\n")

    overall = {'number': 0, 'total': 0, 'shortest':1e6, 'longest': 0}


    for faf in files:
        if args.v:
            sys.stderr.write(f"Counting sequences in {faf}\n")
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

        if args.n:
            continue

        if args.t:
            print("\t".join(map(str, [faf, len(lens), length, lens[0], lens[-1], n50, n75, auN])))
        else:
            print(f"""
File name: {faf}
Number of sequences: {len(lens):,}
Total length: {length:,}
Shortest: {lens[0]:,}
Longest: {lens[-1]:,}
N50: {n50:,}
N75: {n75:,}
auN: {int(auN):,}  """
            )

