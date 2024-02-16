"""
Given an arbitrary CSV file with a column of tax IDs, add the taxonomy.

Note that we throw an error if not all the rows have the same number of columns.
"""

import argparse
import gzip
import sys
from taxon import get_taxonomy_db, get_taxonomy


want = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

dbcon = get_taxonomy_db()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Append taxonomy to a tsv file. We add it as the last column")
    parser.add_argument('-f', help='tsv file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-c', help='taxonomy ID column', required=True, type=int)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    numcols = -1
    
    with gzip.open(args.f, 'rt') as infile, open(args.o, 'w') as out:
        for i, l in enumerate(infile):
            p = l.strip().split("\t")
            if numcols < 0:
                numcols = len(p)
            if len(p) != numcols:
                print("ERROR: There are a different number of columns in this file.", file=sys.stderr)
                print(f"The first line had {numcols} columns. Line {i} has {len(p)} columns", file=sys.stderr)
                sys.exit(1)
            results = ['-', '-', '-', '-', '-', '-', '-', '-', '-']
            t, n = get_taxonomy(p[args.c], dbcon)
            # while t.parent != 1 and t.taxid != 1:
            while t.taxid != 1:
                if t.rank in want and n.scientific_name:
                    results[want.index(t.rank)] = want[want.index(t.rank)][0] + "__" + n.scientific_name
                t, n = get_taxonomy(t.parent, dbcon)
            print("\t".join(map(str, p + results)), file=out)


