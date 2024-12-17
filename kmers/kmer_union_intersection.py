"""
Count the kmers in two genbank files and report the number of k-mers in the
union and intersection of them.
"""

import os
import sys
import argparse

__author__ = 'Rob Edwards'

from Bio import SeqIO

def count_kmers(genbank_file):
    # read the genbank file
    gbk = SeqIO.parse(genbank_file, 'genbank')
    count = {}
    for record in gbk:
        posn = 0
        while posn < len(record.seq) - args.s:
            count[str(record.seq[posn:posn + args.s])] = count.get(record.seq[posn:posn + args.s], 0) + 1
            posn += 1
    return count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='first genbank file', required=True)
    parser.add_argument('-g', help='second genbank file', required=True)
    parser.add_argument('-s', help='kmer size (default = 11)', type=int, default=11)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.v:
        print(f"Counting kmers in {args.f}", file=sys.stderr)
    count_f = count_kmers(args.f)
    set_f = set(count_f.keys())

    if args.v:
        print(f"Counting kmers in {args.g}", file=sys.stderr)
    count_g = count_kmers(args.g)
    set_g = set(count_g.keys())

    intersection = set_f.intersection(set_g)
    union = set_f.union(set_g)

    """
    print(f"The intersection between {args.f} and {args.g} is {len(intersection)}")
    print(f"The union between {args.f} and {args.g} is {len(union)}")
    """
    print(f"{args.f}\t{args.g}\t{args.s}\t{len(intersection)/len(union):.4f}%")