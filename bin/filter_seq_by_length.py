"""
Filter a fasta file for sequences longer than a specified length
"""

import os
import sys
import argparse
from roblib import sequences

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter a fasta file for sequences longer than a specified length')
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-l', help='Minimum length to filter on (seq >= this number)', type=int)
    args = parser.parse_args()

    for seqid, seq in sequences.stream_fasta(args.f):
        if len(seq) >= args.l:
            print(">{}\n{}".format(seqid, seq))