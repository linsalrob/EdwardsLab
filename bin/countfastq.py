#!/usr/bin/env python3

"""
Count the characters in a fasta file. We summarize the longest and shortest reads
and the N50 of the data set.
"""

import os
import sys
import argparse
from roblib import bcolors, stream_fastq

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', nargs='+', help='fastq file')
    parser.add_argument('-d', nargs='+', help='directory of fastq files')
    parser.add_argument('-t', help='tab separated summmary of name, total len, shortest, longest, n50, n75', action="store_true")
    parser.add_argument('-s', help='(deprecated). Same as -t', action="store_true")
    args = parser.parse_args()

    if not args.f and not args.d:
        sys.stderr.write(f"{bcolors.RED}FATAL: Please specify either -f or -d{bcolors.ENDC}\n")
        sys.exit(1)

    if args.f:
        files = args.f
    else:
        files = []

    if args.d:
        for subdir in args.d:
            for f in os.listdir(subdir):
                files.append(os.path.join(subdir, f))

    overall = {'number': 0, 'total': 0, 'shortest':1e6, 'longest': 0}
    for faf in files:
        if not os.path.exists(faf):
            sys.stderr.write(f"{bcolors.RED}FATAL: {faf} not found{bcolors.ENDC}\n")
            sys.exit(1)

        lens = []
        for (sid, label, seq, qual) in stream_fastq(faf):
            lens.append(len(seq))
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

        if args.s or args.t:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(faf, len(lens), length, lens[0], lens[-1], n50, n75))
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
        overall['number'] += len(lens)
        overall['total']  += length
        if lens[0] < overall['shortest']:
            overall['shortest'] = lens[0]
        if lens[-1] > overall['longest']:
            overall['longest'] = lens[-1]


print(f"""

OVERALL SUMMARY
Number of sequences: {overall['number']:,}
Total length: {overall['total']:,}
Shortest: {overall['shortest']:,}
Longest: {overall['longest']:,}
"""
      )



