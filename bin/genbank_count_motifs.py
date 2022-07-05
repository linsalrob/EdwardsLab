"""
Count the occurrence of motifs in a genbank file
"""

import os
import sys
import argparse
from roblib import genbank_seqio, rc


__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count the motif and its occurrence in a sequence')
    parser.add_argument('-f', help='input genbank file', required=True)
    parser.add_argument('-m', help='motif to look for', type=str)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    motif = args.m.upper()

    for seq in genbank_seqio(args.f, args.v):
        dna = seq.seq.upper()
        count = dna.count(motif)
        count += dna.count(rc(motif))
        print(f"{args.f}\t{seq.id}\t{count}")