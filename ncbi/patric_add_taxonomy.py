"""
Add the taxonomy to the patric metadata file
"""

import os
import sys
import argparse
import taxon

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Append taxonomy to the patric metadata file. This adds it at column 67")
    parser.add_argument('-f', help='patric metadata file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-t', help='taxonomy directory (default=/home2/db/taxonomy/current/)',
                        default='/home2/db/taxonomy/current/')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    sys.stderr.write("Reading taxonomy\n")
    taxa = taxon.read_nodes(directory=args.t)
    names, blastname = taxon.read_names(directory=args.t)
    divs = taxon.read_divisions(directory=args.t)

    sys.stderr.write("Read taxonomy\n")
    want = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    with open(args.o, 'w', encoding='utf-8') as out: 
        with open(args.f, 'r', encoding='utf-8') as f:
            for l in f:
                p = l.strip().split("\t")
                while (len(p) <= 68):
                    p.append("")

                if l.startswith("genome_id"):
                    out.write("{}\t{}\n".format(l.strip(), "\t".join(want)))
                    continue

                tid = p[3]

                level = {}
                while tid != '0' and tid != '1' and tid in taxa and taxa[tid].parent != '1':
                    if taxa[tid].rank in want:
                        level[taxa[tid].rank] = names[tid].name
                    tid = taxa[tid].parent


                for w in want:
                    if w in level:
                        p.append(level[w])
                    else:
                        p.append("")

                out.write("\t".join(p))
                out.write("\n")
