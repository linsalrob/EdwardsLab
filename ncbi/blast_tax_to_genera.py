"""
Read a blast file generated with -outfmt '6 std qlen slen staxids sscinames' and generate a list of taxonomies that match.

Our taxonomy has kingdom / phylum / genus / species
"""

import os
import sys
import argparse
import taxon
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read a blast file and create a tuple of [query / kingdom / phylum / genus / species]")
    parser.add_argument('-f', help='blast file(s). Note this must have taxids as column 14. You may specify more than one file', required=True, action='append')
    parser.add_argument('-t', help='taxonomy directory (default=/home2/db/taxonomy/current/)', default='/home2/db/taxonomy/current/')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    want = ['superkingdom', 'phylum', 'genus', 'species']

    sys.stderr.write("Reading taxonomy\n")
    taxa=taxon.read_nodes(directory=args.t)
    names,blastname = taxon.read_names(directory=args.t)
    sys.stderr.write("Read taxonomy\n")

    with open(args.f, 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")

            for tid in p[14].split(";"):
                level = {}
                while tid != '0' and tid != '1' and tid in taxa and taxa[tid].parent != '1':
                    if taxa[tid].rank in want:
                        level[taxa[tid].rank] = names[tid].name
                    tid = taxa[tid].parent

                results = [p[0], tid]
                for w in want:
                    if w in level:
                        results.append(level[w])
                    else:
                        results.append("")

                print("\t".join(results))

