"""
Rename the sequences in a fasta file

If you give an optional -r the sequences will be renamed with that, otherwise with the file name
"""

import os
import sys
import argparse
from roblib import stream_fasta

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rename the sequences in a fasta file')
    parser.add_argument('-f', help='input fasta file', required=True)
    parser.add_argument('-r', help='string to rename the seqids to. Default=use the file name')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.r:
        ren = args.r
    else:
        ren = args.f

    counter = 0
    for seqid, seq in stream_fasta(args.f):
        counter += 1
        print(f">{ren}_{counter}\n{seq}")
