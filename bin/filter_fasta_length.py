"""
Stream a fasta file and print it out
"""

import os
import sys
import argparse
from roblib import sequences
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="stream the contents of a fasta file")
    parser.add_argument('-f', help='file to stream', required=True)
    parser.add_argument('-m', help='minimum sequence length', type=int, default=1000)
    args = parser.parse_args()

    for (seqid, seq) in sequences.stream_fasta(args.f):
        if len(seq) > args.m:
            print(">{}\n{}".format(seqid, seq))
