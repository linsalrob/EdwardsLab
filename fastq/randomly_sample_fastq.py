"""
Randomly sample a fraction of a fastq file and just create one file

Note that if -p is 100 you will randomize the order of sequences in the file
"""

import os
import sys
import argparse
from random import sample
from roblib import stream_fastq, message

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Randomly sample a single fastq file")
    parser.add_argument('-f', help='fastq file to sample', required=True)
    parser.add_argument('-o', help='output file name', required=True)
    parser.add_argument('-p', help='percent of the file to sample', type=int)
    parser.add_argument('-s', help='maxium size to write (Human readable format ok: K, M, G. eg 100M)')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not args.s and not args.p:
        message("Either -s or -p must be provided", "RED")
        exit(-1)

    sequences = []
    for seqid, header, seq, qualscores in stream_fastq(args.f):
        sequences.append([header, seq, qualscores])

    n = len(sequences)
    maxsize = 1e20

    if args.s:
        size = args.s
        # convert the units to 1000's
        size_name = ["B", "K", "M", "G", "T", "P", "E", "Z", "Y"]
        try:
            if size[-1] in size_name:
                unit = size[-1]
                maxsize = (1024 ** size_name.index(unit)) * int(size[:-1])
                message(f"Limiting to {maxsize}", "GREEN")
            else:
                maxsize = int(size)
        except ValueError as e:
            message(f"Can't parse size from {size}", "RED", stderr=True)
            exit(-1)

    if args.p:
        n = int(args.p/100 * len(sequences))

    if args.v:
        message(f"There are {len(sequences)} sequences. So we will sample {n} elements", "GREEN")

    processed_size = 0
    with open(args.o, 'w') as out:
        for s in sample(sequences, n):
            out.write(f"@{s[0]}\n{s[1]}\n+\n{s[2]}\n")
            processed_size += len(s[1])
            if processed_size > maxsize:
                break

