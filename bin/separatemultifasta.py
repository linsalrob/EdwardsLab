#!/usr/bin/env python3

"""
Separate a multi-fasta file into separate files in a directory
"""

import os
import sys
import argparse
from roblib import stream_fasta, colours

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-d', help='output directory', required=True)
    parser.add_argument('-n', help='number of sequences per file (default=1)', default=1, type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.d):
        os.makedirs(args.d, exist_ok=True)


    if args.n == 1:
        for seqid, seq in stream_fasta(args.f, True):
            sname = seqid.split(" ")[0]
            if args.v:
                sys.stderr.write(f"{colours.GREEN}Writing {sname}{colours.ENDC}\n")
            with open(os.path.join(args.d, f"{sname}.fasta"), 'w') as out:
                out.write(f">{seqid}\n{seq}\n")
    else:
        seqsadded = 0
        filecount = 0
        out = open(os.path.join(args.d, f"sequences.{filecount}.fasta"), 'w')
        if args.v:
            sys.stderr.write(f"{colours.GREEN}Writing sequences.{filecount}.fasta{colours.ENDC}\n")
        filecount += 1
        for seqid, seq in stream_fasta(args.f, True):
            out.write(f">{seqid}\n{seq}\n")
            seqsadded += 1
            if seqsadded == args.n:
                seqsadded = 0
                out.close()
                out = open(os.path.join(args.d, f"sequences.{filecount}.fasta"), 'w')
                if args.v:
                    sys.stderr.write(f"{colours.GREEN}Writing sequences.{filecount}.fasta{colours.ENDC}\n")
                filecount += 1
        out.close()
