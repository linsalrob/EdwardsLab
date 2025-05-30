"""
Find runs that match crAssphage with a minimum number of runs, a minimum alignment length, and a minimum percent
identity
"""

import os
import sys
import argparse

import pysam

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count hits in a bam file")
    parser.add_argument('-d', help='Directory of bam files', required=True)
    parser.add_argument('-l', help='Minimum alignment length (default = 30nt)', default=30, type=int)
    parser.add_argument('-c', help='Minimum coverage as a fraction of read length (default = 0.9)', default=0.9, type=float)
    parser.add_argument('-v', help='Verbose output', action='store_true')
    args = parser.parse_args()

    for f in os.listdir(args.d):
        count = {}
        found = False
        if not f.endswith('.bam'):
            if args.v:
                sys.stderr.write(f"Skipped {f}\n")
            continue
        if args.v:
            sys.stderr.write(f"Reading {f}\n")

        bam = pysam.AlignmentFile(os.path.join(args.d, f), 'rb')

        for p in bam.pileup():
            for pilups in p.pileups:
                if pilups.alignment.query_length < args.l:
                    if args.v:
                        sys.stderr.write(f"Query length too short ({pilups.alignment.query_length}) in {f}\n")
                    continue
                alignment_frac = 1.0 * pilups.alignment.query_alignment_length / pilups.alignment.query_length
                if alignment_frac <= args.c:
                    if args.v:
                        sys.stderr.write(f"Alignment fraction too short ({alignment_frac}) in {f}\n")
                    continue
                count[pilups.alignment.reference_name] = count.get(pilups.alignment.reference_name, 0) + 1
                found = True


        if found:
            for q in count:
                print(f"{f}\t{q}\t{count[q]}")
