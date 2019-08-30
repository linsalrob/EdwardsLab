"""
Print the taxonomy associated with an id. If we are given a file, we print each of the items in the file
"""

import os
import sys
import argparse
from taxon import get_taxonomy_db, get_taxonomy

def print_tax(tid, verbose=False):
    """
    Print the taxonomy associated with this id
    :param tid: the taxonomy id
    :param verbose: more output
    :return: nada
    """

    wanted_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
    taxlist = ["", "", "", "", "", "", "", ""]

    # connect to the SQL dataabase
    c = get_taxonomy_db()

    t, n = get_taxonomy(tid, c)
    if not t:
        if verbose:
            sys.stderr.write("No taxonomy for {}\n".format(tid))
        return

    while t.parent != 1 and t.taxid != 1:
        if t.rank in wanted_levels:
            taxlist[wanted_levels.index(t.rank)] = n.scientific_name
        t, n = get_taxonomy(t.parent, c)
    print("\t".join(map(str, [tid]+taxlist)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print taxonomy strings by id")
    parser.add_argument('-f', help='file of taxon ids, one per line')
    parser.add_argument('-t', help='taxonomy ID', action='append') 
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    if not args.t and not args.f:
        sys.stderr.write("One of -t or -f are required.\n")
        sys.exit(-1)

    print("\t".join(['taxid', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']))

    if args.t:
        for t in args.t:
            print_tax(t, args.v)

    if args.f:
        with open(args.f, 'r') as f:
            for l in f:
                print_tax(l.strip(), args.v)
