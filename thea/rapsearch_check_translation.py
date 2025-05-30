"""
Some rapsearch searches have not been translated properly.

This just counts the alphabet in the query: string of the first 2,000 results to see if they have been or not
"""

import os
import sys
import argparse
import gzip

__author__ = 'Rob Edwards'

def read_rapsearch_aln(rpfile, linestocount, verbose=False):
    """
    Read the rapsearch file and just return the query lines
    :param rpfile: rapsearch file to read
    :param linestocount: how many query lines to count
    :param verbose: more output
    :return: a query line
    """
    
    if rpfile.endswith('.gz'):
        qin = gzip.open(rpfile, 'rt')
    else:
        qin = open(rpfile, 'r')

    lines = []
    for l in qin:
        if l.startswith('Query:'):
            lines.append(l)
            if len(lines) > linestocount:
                break

    qin.close()

    return lines

def count_bases(rpfile, linestocount, verbose=False):
    """
    Count the alphabet in the query lines of a rapsearch file
    :param rpfile: the rapsearch aln file
    :param verbose: more output
    :return: a set of the letters
    """

    counts = set()
    counter = 0
    lines = read_rapsearch_aln(rpfile, linestocount, verbose)
    for l in lines:
        p = l.split()
        if len(p) > 4:
            sys.stderr.write("Split wasn't magic on |{}|\n".format(l))
            continue
        counts.update([x.lower() for x in p[2]])

    return counts


def is_dna(counts, rpfile, verbose=False):
    """
    Return true or false for whether a dict contains dna or protein
    :param counts: the set of counts
    :param rpfile: the rapsearch file name (just in case its ambiguous)
    :param verbose: more output
    :return: true or false
    """

    # pure DNA
    for bp in 'acgtn-':
        counts.discard(bp)

    if not counts:
        return "DNA"

    # letters that are ONLY amino acids: E, F, I, L, P, Q

    for aa in 'efilpq':
        if aa in counts:
            return "PROTEIN"

    if verbose:
        sys.stderr.write("AMBIGUOUS. Not sure what to report for {}\n".format(rpfile))

    return "AMBIGUOUS"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='rapsearch output file', required=True)
    parser.add_argument('-l', help='number of lines to search (default = 2000)', default=2000, type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    counts= count_bases(args.f, args.l, args.v)
    print("{}\t{}".format(args.f, is_dna(counts, args.f, args.v)))
