"""
Split the R1 and R2 reads in an interleaved fasta file
"""

import os
import sys
import argparse
from roblib import stream_fasta

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
    parser.add_argument('-1', '--read1', help='R1 output file', required=True)
    parser.add_argument('-2', '--read2', help='R2 output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    with open(args.read1, 'w') as r1out, open(args.read2, 'w') as r2out:
        for seqid, seq in stream_fasta(args.f, True):
            r = seqid.split(" ")[0][-1]
            if r == '1':
                print(f">{seqid}\n{seq}", file=r1out)
            elif r == '2':
                print(f">{seqid}\n{seq}", file=r2out)
            else:
                print(f"Can't parse 1/2 from {seqid}", file=sys.stderr)