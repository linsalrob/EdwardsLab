"""
Split a paired fastqfile into a random set of sequences
"""

import os
import sys
import argparse

from roblib import bcolors
from random import randint
from roblib import bcolors, stream_paired_fastq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot a heatmap")
    parser.add_argument('-l', help='fastq R1 file', required=True)
    parser.add_argument('-r', help='fastq R2 file', required=True)
    parser.add_argument('-o', help='output file stem. Numbers will be appended to this', required=True)
    parser.add_argument('-d', help='output directory', required=True)
    parser.add_argument('-n', help='number of files to split into', required=True, type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.d):
        os.mkdir(args.d)

    files = {}
    for i in range(1, args.n+1):
        if args.v:
            sys.stderr.write(f"{bcolors.GREEN}FILE: {i}{bcolors.ENDC}\n")
        files[i] = [
            open(os.path.join(args.d, args.o + ".R1." + str(i)) + ".fastq", 'w'),
            open(os.path.join(args.d, args.o + ".R2." + str(i)) + ".fastq", 'w')
        ]

    for seqid, header1, seq1, qualscores1, header2, seq2, qualscores2 in stream_paired_fastq(args.l, args.r):
        outint = randint(1, args.n)
        if args.v:
            sys.stderr.write(f"{bcolors.PINK}FILE: {outint}{bcolors.ENDC}\n")
        files[outint][0].write(f"@{header1}\n{seq1}\n+\n{qualscores1}\n")
        files[outint][1].write(f"@{header2}\n{seq2}\n+\n{qualscores2}\n")

    for i in range(1, args.n + 1):
        files[i][0].close()
        files[i][1].close()

