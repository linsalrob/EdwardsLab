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
    filecount = 0
    for dirpath, dirname, filename in os.walk(args.d):
        for paf in filename:
            if args.v:
                sys.stderr.write(f"Processing {os.path.join(dirpath, paf)}\n")
            filecount += 1
            with gzip.open(os.path.join(dirpath, paf), 'rt') as fin:
                for l in fin:
                    p = l.strip().split("\t")
                    if p[5] not in counts:
                        counts[p[5]] = [0, 0]
                    counts[p[5]][0] += int(p[9])
                    counts[p[5]][1] += int(p[10])

    sys.stderr.write(f"Processed {filecount} PAF files\n")
    with open(args.o, 'w') as out:
        for contig in counts:
            out.write(f"{contig}\t{counts[contig][0]}\t{counts[contig][1]}\t{counts[contig][0]/counts[contig][1]}\n")

