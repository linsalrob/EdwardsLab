"""
Reservior sample a fasta file

This is a single pass (c.f. subsample_fasta.py which uses 2 passes) as the reservior only needs that many!
"""
import gzip
import os
import sys
import argparse
from roblib import stream_fasta, bcolors, message
import random

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', '--fasta', help='fasta file', required=True)
    parser.add_argument('-r', '--reverse', help='r2 read file')
    parser.add_argument('-o', '--output', help='output fasta file', required=True)
    parser.add_argument('-p', '--reverse_output', help='')
    parser.add_argument('-s', '--sample', help='how many sequences to subsample (default=100,000)',
                        type=int, default=100000)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.reverse and args.reverse_output is None:
        parser.error("The --reverse_output option is required if you include a --reverse option")

    # the reservoirs
    ids = []
    seqs = []
    if args.verbose:
        print(f"{bcolors.BLUE}Parsing {args.fasta}{bcolors.ENDC}", file=sys.stderr)
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

    if args.reverse:
        finalids = set(ids)

        ron = 0
        if args.verbose:
            print(f"{bcolors.BLUE}Parsing {args.reverse}{bcolors.ENDC}", file=sys.stderr)
        with opener(args.reverse_output, 'wt') as out:
            for seqid, seq in stream_fasta(args.reverse, whole_id=True):
                if seqid in finalids:
                    print(f">{seqid}\n{seq}", file=out)
                    ron += 1
        if ron != len(finalids):
            message(f"WARNING: We printed {len(finalids)} in the forward file but only {ron} ids in the reverse", "RED")




