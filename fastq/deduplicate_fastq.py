"""
Remove duplicate reads from a fastq file and check their sequences are the same
"""

import os
import sys
import argparse
from roblib import stream_fastq, bcolors

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='fastq file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    with open(args.o, 'w') as out:
        seen = {}
        for sid, header, seq, qual in stream_fastq(args.f):
            if sid in seen:
                if seen[sid] != seq:
                    print(f"{bcolors.RED}ERROR: {sid} has different sequences: =======\n{seen[sid]}\n\n{seq}\n======={bcolors.ENDC}", file=sys.stderr)
                    continue
                if args.v:
                    sys.stderr.write(f"Duplicate sequence {sid} ignored\n")
            else:
                seen[sid] = seq
                out.write(f"@{header}\n{seq}\n+\n{qual}\n")
