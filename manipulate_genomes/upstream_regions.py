"""
Extract the region from the beginning of the genome (position 1) to the start of the first gene.

This is really desgined for picobirnaviruses
"""

import os
import sys
import argparse
from roblib import genbank_seqio

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract upstream sequences')
    parser.add_argument('-f', help='input genbank file', required=True)
    parser.add_argument('-x', help='maximum upstream sequence', type=int, default=100000000)
    parser.add_argument('-m', help='minimum upstream sequence', type=int, default=0)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    for seq in genbank_seqio(args.f, args.v):
        firststart = 100000000
        strand = 0
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            (start, stop, strand) = (feat.location.start.position, feat.location.end.position, feat.strand)
            if start < firststart:
                firststart = start
                strand = 0
            if stop < firststart:
                firststart = stop
                strand = -1
        if strand < 0:
            sys.stderr.write(f"First position is on -ve strand for {seq.id} so ATG not present\n")
        firststart += 3
        fragseq = seq.seq[0:firststart]
        if len(fragseq) > args.x:
            fragseq = fragseq[-args.x:]
        if len(fragseq) > args.m:
            print(f">{seq.id} [{seq.description}]\n{fragseq}")