"""
Read a single protein fasta files and calculate md5sums for proteins

NOTE: This is faster (not really fast), but only allows you to run one fasta
"""

import os
import sys
import argparse
import hashlib
from roblib import stream_fasta, bcolors


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='fasta file')
    parser.add_argument('-i', help='id map file to write', required=True)
    parser.add_argument('-o', help='output file', required=True)
    args = parser.parse_args()

    md5seen = set()
    idseen = set()

    with open(args.o, 'w') as out, open(args.i, 'w') as idout:
        for seqid, seq in stream_fasta(args.f):
            if seqid in idseen:
                continue
            idseen.add(seqid)
            md5 = hashlib.md5(seq.upper().encode('utf-8')).hexdigest()
            idout.write(f"{md5}\t{seqid}\n")
            if md5 not in md5seen:
                out.write(f">{md5}\n{seq}\n")
            md5seen.add(md5)

