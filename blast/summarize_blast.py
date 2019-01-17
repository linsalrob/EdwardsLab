"""
Summarize hits in a blast file based on the sequences present
"""

import os
import sys
import argparse
from roblib import stream_fasta, stream_blast_results

def seq_lengths(fafile, verbose=False):
    """
    Read the sequence length from a fasta file
    :param fafile: the fasta file to read
    :param verbose: more output
    :return: a dict of sequence id and length
    """

    length = {}
    for i,s in stream_fasta(fafile):
        length[i] = len(s)

    return length


def summarize_blast(fafile, blfile, verbose=False):
    """
    Summarize blast hits
    :param fafile: the query fasta file
    :param blfile: the blast output file
    :param verbose: more output
    :return:
    """

    seqlens = seq_lengths(fafile, verbose)

    for b in stream_blast_results(blfile, verbose=verbose):
        "blech blech blech"





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', help='blast output file', required=True)
    parser.add_argument('-f', help='fasta query file', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

