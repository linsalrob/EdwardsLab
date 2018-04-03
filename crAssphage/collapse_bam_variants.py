"""
Collapse variants from a bamfile and try and make as few variants as possible
"""

import os
import sys
import argparse
from roblib import read_fasta

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Collapse all variants from a bam file and make as few options as possible")
    parser.add_argument('-b', help='bam file', required=True)
    parser.add_argument('-r', help='reference fasta sequence', required=True)
    parser.add_argument('-s', help='start of the region to look at (default = whole sequence)', default=0, type=int)
    parser.add_argument('-e', help='end of sequence to look at (default = whole sequence)', type=int)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    fa = read_fasta(args.r)

    if not args.e:
        args.e = max([len(fa[f]) for f in fa])



