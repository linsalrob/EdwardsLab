"""
Get all taxids from the database and make a sorted list of tid/rank
"""

import os
import sys
import argparse
from taxon import get_taxonomy_db, get_taxonomy, all_ids

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

def find_rank(tid, trank, tdb, verbose=False):
    """
    Find the value for trank starting at tid.
    :param tid: The taxonomy ID
    :param trank: The taxonomic rank to return
    :param tdb: The taxonomy database
    :param verbose: More output
    :return: the taxonomic rank for tid or root if it is not found
    """
    global id2rank
    if tid in id2rank:
        return id2rank[tid]

    seenids = set()

    t,n = get_taxonomy(tid, tdb)
    while t.parent != 1 and t.taxid != 1 and t.rank != trank and t.taxid not in id2rank:
        seenids.add(t.taxid)
        t, n = get_taxonomy(t.parent, tdb)
    if t.taxid in id2rank:
        for s in seenids:
            id2rank[s] = id2rank[t.taxid]
        return id2rank[t.taxid]
    if t.rank == trank:
        for s in seenids:
            id2rank[s] = n.scientific_name
        id2rank[tid] = n.scientific_name
        return n.scientific_name
    if verbose:
        sys.stderr.write(f"{colours.PINK}ERROR: No rank for {tid}\n{colours.ENDC}")

    return "root"




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print the taxonomic rank for all ids ")
    parser.add_argument('-t', help='taxnomic rank', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    c = get_taxonomy_db()
    for i in all_ids(c, args.v):
        print(f"{i}\t{find_rank(i, args.t, c, args.v)}")
    