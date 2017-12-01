"""
Create a python taxonomy for Mike Doane based on a tsv file
"""

import os
import sys
import argparse
import taxon

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse a tsv file and add taxonomy')
    parser.add_argument('-f', help='tab seperated values', required=True)
    args = parser.parse_args()

    want = ['superkingdom', 'kingdom' 'phylum', 'class', 'order', 'family', 'genus', 'species']

    sys.stderr.write("Reading databases\n")
    taxa = taxon.read_nodes()
    names, blastname = taxon.read_names()
    sys.stderr.write("Done\n")

    with open(args.f, 'r') as f:
        for l in f:
            if l.startswith("#"):
                print("{}\t".format(l.strip()) + "\t".join(want))
                continue
            p = l.strip().split("\t")
            m = ["" for w in want]
            i = p[2]
            while taxa[i].parent != '1' and i != '1':
                if taxa[i].rank in want:
                    m[want.index(taxa[i].rank)] = names[i].name
                i = taxa[i].parent
            print("{}\t".format(l.strip()) + "\t".join(m))
