"""
Filter a fasta file on length
"""

import os
import sys
import argparse
from roblib import stream_fasta

def length_filter(f, l, verbose=False):
    """
    Filter a fasta file based on the minimum length, l
    :param f: fasta file
    :param l: minimum sequene length
    :param verbose: more output
    :return:
    """

    for seqid, seq in stream_fasta(f, True):
        if len(seq) < l:
            continue
        print(">{}\n{}".format(seqid, seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter a file based on length")
    parser.add_argument('-f', help='the fasta file to filter', required=True)
    parser.add_argument('-l', help='minimum length (default=1000)', default=1000, type=int)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    length_filter(args.f, args.l, args.v)