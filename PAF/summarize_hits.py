"""
PAF files have the number of hits in col 10 and the alignment len in column 11.
This script finds all paf files in a directory (recursively) and then summarizes all hits

At the moment it does not consider location of the hits!
"""

import os
import sys
import argparse
import gzip

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summarize PAF hits')
    parser.add_argument('-d', help='directory with paf files (will be recursed)', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    counts = {}
    for dirpath, dirname, filename in os.walk(args.d):
        for paf in filename:
            if args.v:
                sys.stderr.write(f"Processing {os.path.join(dirpath, paf)}")
            with gzip.open(os.path.join(dirpath, paf), 'rt') as fin:
                for l in fin:
                    p = l.strip().split("\t")
                    if p[5] not in counts:
                        counts[p[5]] = [0, 0]
                    counts[p[5]][0] += int(p[9])
                    counts[p[5]][1] += int(p[10])

    for contig in counts:
        print(f"{contig}\t{counts[0]}\t{counts[1]}\t{counts[contig][0]/counts[contig][1]}")

