"""
Take any line with a patric ID and prepend taxonomy
"""

import os
import sys
import argparse
from taxon import connect_to_db, get_taxonomy, get_taxid_for_name
import re

__author__ = 'Rob Edwards'

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

taxonomy_str = {}

def resolve_taxonomy(tid, conn, verbose=False, nocolor=False):
    """
    Convert the taxonomy id to a tab separated string
    :param tid: the taxonomy object
    :param conn: the database connection
    :param verbose: more output
    :param nocolor: no color ouput
    :return: a string representing the taxonomy
    """

    global taxonomy_str
    if tid in taxonomy_str:
        return taxonomy_str[tid]

    wanted_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    rnk = ['', '', '', '', '', '', '']
    t, n = get_taxonomy(tid, conn)
    while t.parent != 1 and t.taxid != 1:
        if t.rank in wanted_levels:
            rnk[wanted_levels.index(t.rank)] = n.scientific_name
        t, n = get_taxonomy(t.parent, conn)
    taxonomy_str[tid] = "\t".join(map(str, rnk))
    return taxonomy_str[tid]

def add_tax(infile, verbose=False):
    """
    Add the taxonomy to any file with a PATRIC id
    :param infile: the strain results file from Greg
    :param verbose: more output
    :return:
    """

    db = connect_to_db("/data/ncbi/taxonomy.sqlite3")

    with open(infile, 'r') as f:
        for l in f:
            if 'PATRIC' in l:
                m=re.search('PATRIC\|(\d+)\.\d+', l)
                tid = m.groups(0)[0]
                sys.stderr.write(f"{bcolors.OKGREEN}Resolving: {tid}{bcolors.ENDC}\n")
                tax = resolve_taxonomy(tid, db, args.v)
                sys.stdout.write(f"{tax}\t{l}")
            else:
                sys.stdout.write(f"{l}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='file to parse', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    add_tax(args.f, args.v)