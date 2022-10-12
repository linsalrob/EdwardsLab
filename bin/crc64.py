"""
Calculate the crc64 checksum of a fasta sequence.
"""

import os
import sys
import argparse
from crc64iso.crc64iso import crc64
from roblib import stream_fasta

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-o', help='output file (default -)', default=sys.stdout)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    for seqid, seq in stream_fasta(args.f):
        print(f"{seqid}\t{crc64(seq)}")