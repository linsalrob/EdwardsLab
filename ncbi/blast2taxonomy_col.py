"""

"""

import os
import sys
import argparse

from roblib import bcolors
from taxon import get_taxonomy_db, get_taxonomy, get_taxid_for_name, taxonomy_hierarchy

taxa = {}

def taxstring(tid, verbose=False):
    """

    :param tid: taxonomy ID
    :param verbose: more output
    :return: an array of the taxnomy from kingdom -> species
    """
    global taxa
    if tid in taxa:
        return taxa[tid]

    want = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    thistaxa = ['', '', '', '', '', '', '']
    c = get_taxonomy_db()
    for p in taxonomy_hierarchy(tid, verbose=False):
        m,n = get_taxonomy(p, c)
        if m.rank in want:
            thistaxa[want.index(m.rank)] = n.blast_name[0].upper() + n.blast_name[1:]
    taxa[tid] = thistaxa
    return taxa[tid]


def parse_blast(bf, taxcol, verbose=False):
    """

    :param bf: the blast output file
    :param taxcol: the column that contains the taxonomy ID
    :param verbose: more output
    :return:
    """

    lastcol = -1
    with open(bf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if lastcol == -1:
                lastcol = len(p)
            if len(p) != lastcol:
                sys.stderr.write(f"{bcolors.RED}FATAL: Uneven number of columns. We had {lastcol} but now {len(p)}\n")
                sys.exit(-1)

            t = taxstring(p[taxcol], verbose)
            print("\t".join(p+t))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', help='blast output file', required=True)
    parser.add_argument('-c', help='column that has the taxonomy ID', required=True, type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    parse_blast(args.f, args.c, args.v)