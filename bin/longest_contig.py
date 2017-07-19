"""
Print the length of the longest contig for each file in a directory of fasta files.
"""

import os
import sys
import argparse
from roblib import read_fasta

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print the length of the longest contig for each file in a directory of fasta files')
    parser.add_argument('-d', help='Directory of fasta files', required=True)
    args = parser.parse_args()

    for f in os.listdir(args.d):
        fa = read_fasta(os.path.join(args.d, f))
        lengths = [len(fa[x]) for x in fa]
        lengths.sort()
        print("{}\t{}".format(f, lengths[-1]))


