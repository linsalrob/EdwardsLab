"""
Randomly sample a fastq file
"""

import os
import sys
import argparse

from roblib import bcolors
from random import randint
from roblib import bcolors, stream_fastq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Randomly sample a single fastq file")
    parser.add_argument('-f', help='fastq file to sample', required=True)
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
        files[i] = open(os.path.join(args.d, args.o + "." + str(i)) + ".fastq", 'w')

    for seqid, header, seq, qualscores in stream_fastq(args.f):
        outint = randint(1, args.n)
        if args.v:
            sys.stderr.write(f"{bcolors.PINK}FILE: {outint}{bcolors.ENDC}\n")
        files[outint].write(f"@{header}\n{seq}\n+\n{qualscores}\n")

    for i in range(1, args.n + 1):
        files[i].close()

