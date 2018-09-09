"""
Create a table of distributions of sequence lengths from one or more fasta files
"""

import os
import sys
import argparse
from roblib import stream_fasta

__author__ = 'Rob Edwards'


def count_len(fastaf, verbose=False):
    """
    Count the sequence lengths and return a dict of len:count
    :param fastaf: fasta file
    :param verbose: more output
    :return:
    """

    counts = {}
    for seqid, seq in stream_fasta(fastaf):
        counts[len(seq)] = counts.get(len(seq), 0) + 1

    return counts

def print_lens(lengths, verbose=False):
    """
    Print all the lengths
    :param lengths: dict of dicts
    :param verbose: more output
    :return:
    """

    lenset = set()
    for x in lengths:
        lenset.update([k for k in lengths[x].keys()])

    alllens = sorted(list(lenset))
    fafs = lengths.keys()

    print("\t{}".format("\t".join(fafs)))
    for l in alllens:
        sys.stdout.write("{}".format(l))
        for f in fafs:
            sys.stdout.write("\t{}".format(lengths[f].get(l, 0)))
        sys.stdout.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Distributions of sequence lengths from one or more fasta files')
    parser.add_argument('-f', help='fasta file(s). Use multiple -f for multiple files', required=True, action='append')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    lengths = {}
    for faf in args.f:
        lengths[faf] = count_len(faf)
    print_lens(lengths)
