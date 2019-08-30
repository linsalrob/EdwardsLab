"""

"""

import os
import sys
import argparse
import gzip
from roblib import bcolors
from taxon import get_taxonomy_db, get_taxonomy, taxonomy_hierarchy, Error

taxa = {}

def choosename(n, verbose=True):
    if n.blast_name:
        return n.blast_name
    if n.scientific_name:
        return n.scientific_name
    if n.common_name:
        return n.common_name
    if n.equivalent_name:
        return n.equivalent_name
    return None

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
    m, n = get_taxonomy(tid,c)
    thisname = choosename(n, verbose)
    if thisname:
        if m.rank in want:
            thistaxa[want.index(m.rank)] = thisname[0].upper() + thisname[1:]
    for p in taxonomy_hierarchy(tid, verbose=False):
        m,n = get_taxonomy(p, c)
        thisname = choosename(n, verbose)
        if not thisname:
            sys.stderr.write(f"{bcolors.RED}ERROR: No name for {tid}{bcolors.ENDC}\n")
            return
        if m.rank in want:
            thistaxa[want.index(m.rank)] = thisname[0].upper() + thisname[1:]
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
    if bf.endswith('.gz'):
        f = gzip.open(bf, 'rt')
    else:
        f = open(bf, 'r')
    for l in f:
        p = l.strip().split("\t")
        if lastcol == -1:
            lastcol = len(p)
        if len(p) != lastcol:
            sys.stderr.write(f"{bcolors.RED}FATAL: Uneven number of columns. We had {lastcol} but now {len(p)}\n")
            sys.exit(-1)

        t = taxstring(p[taxcol], verbose)
        print("\t".join(p+t))
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', help='blast output file', required=True)
    parser.add_argument('-c', help='column that has the taxonomy ID', required=True, type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    parse_blast(args.f, args.c, args.v)
