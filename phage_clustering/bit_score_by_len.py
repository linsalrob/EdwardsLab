"""

"""

import os
import sys
import argparse
from roblib import bcolors, stream_blast_results
__author__ = 'Rob Edwards'



def bit_scores_len(blastf, verbose=False):
    """
    Generate a dict of self:self bitscores
    """

    for b in stream_blast_results(blastf, verbose):
        if b.query == b.db:
            print(f"{b.query_length}\t{b.bitscore}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-b', help='blast input file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    bit_scores_len(args.b, args.v)
