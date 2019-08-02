"""
Read a list of phyla and print out their kindga
"""

import os
import sys
import argparse
from roblib import bcolors
from taxon import get_taxonomy_db, get_taxonomy, get_taxid_for_name, taxonomy_hierarchy
__author__ = 'Rob Edwards'


def print_kingdom(f, verbose=False):

    c = get_taxonomy_db()
    with open(f, 'r') as fin:
        for l in fin:
            l = l.strip()
            t = get_taxid_for_name(l, c)
            if not t:
                print(f"{l}\tUnknown")
                continue
            for p in taxonomy_hierarchy(t, verbose=False):
                m,n = get_taxonomy(p, c)
                if m.rank == 'superkingdom':
                    nm = n.blast_name[0].upper() + n.blast_name[1:]
                    print(f"{l}\t{nm}")





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    
    print_kingdom(args.f, args.v)

