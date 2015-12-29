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
    args = parser.parse_args()

    for (seqid, seq) in sequences.stream_fasta(args.f):
        print("{}\t{}".format(seqid, seq))