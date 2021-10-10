"""
Given a file with tuples of [md5 sum, patric id and function] e.g.

cbcbf68d89875fdd5bc7148798dc7fc2        fig|1560201.3.peg.473 hypothetical protein

We calculate the LCA for each md5 sum

You can then feed that into mmseqs to make a taxonomy database!

At the moment this is not really done in a memory sensitive way :)
"""

import os
import sys
import argparse
import re

from roblib import message
from taxon import get_taxonomy_db, taxonomy_ids_as_list, get_taxonomy

__author__ = 'Rob Edwards'

taxonomies = {}
ignore_taxa = set()


def lca(tids, verbose=False):
    """
    Calculate the lowest common ancestor for a set of taxonomy IDs
    :param set tids: a set of taxonomy IDs
    :param bool verbose: more output
    :return int: the taxonomy ID of the lowest set
    """

    if len(tids) == 1:
        return tids.pop()

    global taxonomies
    global ignore_taxa

    conn = get_taxonomy_db()

    # just for reference
    wanted_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxas = []
    for tid in tids:
        if tid in ignore_taxa:
            continue
        if tid not in taxonomies:
            t, n = get_taxonomy(tid, conn)
            if not t:
                ignore_taxa.add(tid)
                continue
            taxonomies[tid] = taxonomy_ids_as_list(conn, tid, verbose)
        taxas.append(taxonomies[tid])

    # now we just have to add the appropriate taxa and return
    # this is a two dimensional list
    counts = [set(), set(), set(), set(), set(), set(), set()]
    for t in taxas:
        for i, j in enumerate(t):
            counts[i].add(j)
    for i in range(6, -1, -1):
        if len(counts[i]) == 1:
            # this is the lca!
            return counts[i]
    return 1 # root!



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='file of protein MD5 and IDs', required=True)
    parser.add_argument('-o', help='output file of MD5 and taxonomy ID', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    # read all the md5s and map to a list of the genomes
    md5 = {}
    fid = re.compile(r'^fig\|(\d+)\.\d+')
    with open(args.f, 'r') as fin:
        for l in fin:
            p = l.strip().split("\t")
            m = fid.search(p[1])
            if not m:
                message(f"ERROR: Can't parse fig id from {l}", color="RED")
                continue
            if p[0] not in md5:
                md5[p[0]] = set()
            md5[p[0]].add(int(m.groups()[0]))

    with open(args.o, 'w') as out:
        for m in md5:
            tid = lca(md5[m])
            out.write(f"{m}\t{tid}\n")
