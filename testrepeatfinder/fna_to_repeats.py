"""
Read a fasta file and convert it to a table of repeats
"""

import os
import sys
import argparse

from roblib import bcolors, stream_fasta
import RobRepeatFinder

def get_reps(faf, outf, verbose=False):
    """
    Get the repeats and write them out
    :param faf: fasta file
    :param outf: output file
    :param verbose: more output
    :return:
    """

    with open(outf, 'w') as out:
        for seqid, seq in stream_fasta(faf):
            for r in RobRepeatFinder.repeatFinder(seq, 3):
                out.write(f"{r['first_start']}\t{r['first_end']}\t{r['second_start']}\t{r['second_end']}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot a heatmap")
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    get_reps(args.f, args.o, args.v)