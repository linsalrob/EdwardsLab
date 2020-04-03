"""
Extract a sequence from a fasta file
"""

import os
import sys
import argparse
from roblib import stream_fasta



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-i', help='sequence id, (multiple allowed)', nargs='+')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    for seqid, seq in stream_fasta(args.f):
        if seqid in args.i:
            print(f">{seqid}\n{seq}\n")
