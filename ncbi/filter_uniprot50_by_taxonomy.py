"""
The uniref50 data set is available here: ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/

The fasta lines have the taxid of the most common organism (LCA) and so here we filter by taxonomy

"""

import os
import sys
import argparse
import re

from roblib import colour, stream_fasta
from taxon import get_taxonomy_db, get_taxonomy

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

want = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']


def find_rank(tid, trank, tdb, verbose=False):
    """
    Find the value for trank starting at tid.
    :param tid: The taxonomy ID
    :param trank: The taxonomic rank to return
    :param tdb: The taxonomy database
    :param verbose: More output
    :return: the taxonomic rank for tid or root if it is not found
    """

    t,n = get_taxonomy(tid, tdb)
    while t.parent != 1 and t.taxid != 1 and t.rank != trank:
        t, n = get_taxonomy(t.parent, tdb)
    if t.rank == trank:
        return n
    if verbose:
        sys.stderr.write(f"{colour.PINK}ERROR: No rank for {tid}\n{colour.ENDC}")

    return "root"

def fasta2tax(faf, outdir, trank, tdb, verbose=False):
    """
    Split the fasta file by taxonomy
    :param faf: fasta file to split
    :param outdir: output directory to write to
    :param trank: taxonomic rank to choose
    :param tdb: The taxonomy database
    :param verbose: more output
    :return:
    """

    s = re.compile('TaxID=(\d+)')

    fhs = {}
    for seqid, seq in stream_fasta(faf):
        m = s.search(seqid)
        if m:
            tid = m.groups()[0]
            rnk = find_rank(tid, trank, tdb, verbose)
            if rnk not in fhs:
                fhs[rnk] = open(os.path.join(outdir, rnk + ".fasta"), 'w')
            fhs[rnk].write(f">{seqid}\n{seq}\n")
        else:
            sys.stderr.write(f"{colour.RED}ERROR: No taxonomy in {seqid}{colour.ENDC}\n")

    for fh in fhs:
        fhs[fh].close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='uniref file to filter', required=True)
    parser.add_argument('-d', help='output directory', required=True)
    parser.add_argument('-t', required=True,
                        help='Taxonomic rank should probably be one of: ' . ';'.join(want))
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    tdb = get_taxonomy_db()

    fasta2tax(args.f, args.d, args.t, tdb, args.v)