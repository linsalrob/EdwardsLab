"""
Given an accession, get the taxonomy
"""

import os
import sys
import argparse
from roblib import bcolors
from taxon import get_taxonomy_db, acc_to_taxonomy
__author__ = 'Rob Edwards'


def get_taxonomy(acc, protein=False, verbose=False):
    c = get_taxonomy_db()
    return acc_to_taxonomy(acc, c, protein=protein, verbose=verbose)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert a Nucleotide or Protein accession to taxonomy')
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('-a', help='one accession e.g. NC_048037')
    grp.add_argument('-f', help='a file of ids, one per line')
    parser.add_argument('-p', help='accession is for a protein', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.a:
        p,t,n = get_taxonomy(args.a, args.p, args.v)
        print(f"p: {p}\nt: {str(t)}\nn: {n}")
    elif args.f:
        with open(args.f, 'r') as fin:
            for l in fin:
                l = l.strip()
                try:
                    p,t,n = get_taxonomy(l, args.p, args.v)
                except:
                    print(f"{l}\t")
                    print(f"ERROR: There was an error getting the taxonomy ID for {l}", file=sys.stderr)
                    continue

                if t:
                    print(f"{l}\t{t.taxid}")
                else:
                    print(f"{l}\t")
    else:
        print("ERROR: Please choose -a or -f. Use -h for more information", file=sys.stderr)

    


