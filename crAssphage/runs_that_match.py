"""
Find runs that match crAssphage with a minimum number of runs, a minimum alignment length, and a minimum percent
identity
"""

import os
import sys
import argparse

import pysam

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find bam files that match crAssphage")
    parser.add_argument('-d', help='Directory of bam files', required=True)
    parser.add_argument('-l', help='Minimum alignment length (default = 30nt)', default=30, type=int)
    parser.add_argument('-c', help='Minimum coverage as a fraction of read length (default = 0.9)', default=0.9, type=float)
    parser.add_argument('-n', help='Minimum number of reads that meet these criteria (default = 2)', default=2, type=int)
    parser.add_argument('-v', help='Verbose output', action='store_true')
    args = parser.parse_args()

    for f in os.listdir(args.d):
        count = 0
        if not f.endswith('.bam'):
            if args.v:
                sys.stderr.write("SKipped {}\n".format(f))
            continue
        bam = pysam.AlignmentFile(os.path.join(args.d, f), 'rb')

        for p in bam.pileup(truncate=True):
            for pilups in p.pileups:
                if pilups.alignment.query_length < args.l:
                    continue
                alignment_frac = 1.0 * pilups.alignment.query_alignment_length / pilups.alignment.query_length
                if alignment_frac <= args.c:
                    continue
                count += 1
                if count >= args.n:
                    break
            if count >= args.n:
                break
        if count >= args.n:
            print(f)
