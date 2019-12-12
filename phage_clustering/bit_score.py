"""

"""

import os
import sys
import argparse
from roblib import bcolors, stream_blast_results
__author__ = 'Rob Edwards'

def self_bit_scores(blastf, verbose=False):
    """
    Generate a dict of self:self bitscores
    """

    ss = {}
    for b in stream_blast_results(blastf, verbose):
        if b.query == b.db:
            if b.query in ss and ss[b.query] > b.bitscore:
                continue
            ss[b.query] = b.bitscore
    return ss


def pairwise_bit_scores(blastf, ss, verbose=False):
    """
    Make a pairwise average bit score that is 
    the bitscore / average of two proteins self/self bit
    score
    ;param blastf: the blastfile
    :param ss: the self-self bitscores 
    :param verbose: more output
    :return a dict of all vs. all normalized bit scores
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Creating scores{bcolors.ENDC}\n")

    pb = {}
    for b in stream_blast_results(blastf, verbose):
        if b.query not in pb:
            pb[b.query] = {}
        if b.db not in pb:
            pb[b.db] = {}

        # we normalize by the bitscore of the two proteins if we can!
        if b.query in ss and b.db in ss:
            nb = b.bitscore / ((ss[b.query] + ss[b.db])/2)
        else:
            # if we can't do that, we cheat and normalize 
            # the bit score by twice
            # the average length of the proteins
            # i.e. the sum of the lengths
            nb = b.bitscore / (b.query_length + b.subject_length + 3.3)

        if b.query in pb[b.db] and pb[b.db][b.query] > nb:
            continue
        pb[b.db][b.query] = pb[b.db][b.query] = nb
        print(f"{b.query}\t{b.db}\t{b.query_length}\t{b.subject_length}\t{b.bitscore}\t{nb}")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-b', help='blast input file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    ss = self_bit_scores(args.b, args.v)
    pairwise_bit_scores(args.b, ss, args.v)
