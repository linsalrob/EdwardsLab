"""
Print the IDs and lengths of sequences in a fasta file
"""

import os
import sys
import argparse
from roblib import stream_fasta

parser = argparse.ArgumentParser(description="Print the lengths of sequences in a fasta file")
parser.add_argument('-f', help='fasta file', required=True)
parser.add_argument('-w', help='whole sequence ID. Default is to use ID upto whitespace', action="store_true", default=False)
parser.add_argument('-v', help='verbose output', action="store_true")
args = parser.parse_args()

for seqid, seq in stream_fasta(args.f, args.w):
    print("{}\t{}").format(seqid, len(seq))