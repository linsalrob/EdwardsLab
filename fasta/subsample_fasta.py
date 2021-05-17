"""
Subsample a fasta file to generate a number of sequences.

We iterate the file twice to count the number of sequences.
"""

import os
import sys
import argparse
from roblib import stream_fasta, message
from random import shuffle

def get_sequence_ids(inputfile, verbose=False):
    """
    count the sequences in the file
    :param inputfile: the input file 
    :param verbose:  more output
    :return: an array of sequence IDs
    """
    seqids = []
    for seqid, seq in stream_fasta(inputfile, True):
        seqids.append(seqid)

    shuffle(seqids)

    if verbose:
        message(f"There are {len(seqids)} sequences in inputfile", "GREEN")

    return seqids

def subsample(inputfile, n, seqids, verbose=False):
    """
    subsample n sequences from inputfile with counter total sequences
    :param inputfile: fasta file
    :param n: number of sequences
    :param seqids: the array of sequence ids
    :param verbose: more output
    """

    towrite = set(seqids[0:n])

    written = 0
    for seqid, seq in stream_fasta(inputfile, True):
        if seqid in towrite:
            written += 1
            print(f">{seqid}\n{seq}")
    if verbose:
        message(f"Wrote {written} sequences", "GREEN")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-n', help='number of sequences to write', required=True, type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    seqids = get_sequence_ids(args.f, args.v)
    subsample(args.f, args.n, seqids, args.v)
