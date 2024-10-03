"""
Filter out metagenomes

I have created a set of files based on the alignment reports that include only the Mycobacteria reads. See /home/edwa0468/Projects/FAME/SAGCQA0545_combined/Mycobacteria

See https://edwards.flinders.edu.au/percent-similarity-at-different-taxonomic-levels/ to decide on percent ID

Filter those to find a set of patients with "reliable" metagenomes
"""
import gzip
import os
import sys
import argparse

__author__ = 'Rob Edwards'

from pathlib import Path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', '--directory', help='base directory with the alignment files', required=True)
    parser.add_argument('-p', '--percent', help='minimum percent identity for a match (default: 90%)',
                        type=float, default=0.9)
    parser.add_argument('-n', '--number', help='minimum number of proteins to match in this metagenome (default: 2)',
                        type=int, default=2)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.percent > 1:
        args.percent /= 100

    # find the files

    for path in Path(args.directory).rglob('*.tophit_aln.gz'):
        metagenome = os.path.split(path)[-2]
        count = {}
        with gzip.open(path, 'rt') as f:
            for l in f:
                p = l.strip().split("\t")
                if float(p[2]) < args.percent:
                    continue
                count[p[0]] = count.get(p[0], 0) + 1
        if len(count.keys()) > args.number:
            print("\t".join(map(str, [metagenome, path, len(count.keys())])))
