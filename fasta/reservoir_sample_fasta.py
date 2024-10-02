"""
Reservior sample a fasta file

This is a single pass (c.f. subsample_fasta.py which uses 2 passes) as the reservior only needs that many!
"""
import gzip
import os
import sys
import argparse
from roblib import stream_fasta, bcolors
import random

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', '--fasta', help='fasta file', required=True)
    parser.add_argument('-o', '--output', help='output fasta file', required=True)
    parser.add_argument('-s', '--sample', help='how many sequences to subsample (default=100,000)',
                        type=int, default=100000)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    # the reservoirs
    ids = []
    seqs = []
    stream = stream_fasta(args.fasta, whole_id=True)
    for i in range(args.sample):
        seqid, seq = next(stream)
        ids.append(seqid)
        seqs.append(seq)

    n = args.sample
    for seqid, seq in stream:
        j = random.randint(0, n)
        if j < args.sample:
            ids[j] = seqid
            seqs[j] = seq
        n += 1

    opener = open
    if args.output.endswith('.gz'):
        opener = gzip.open
    with opener(args.output, 'wt') as out:
        for i in range(args.sample):
            print(f">{ids[i]}\n{seqs[i]}", file=out)


