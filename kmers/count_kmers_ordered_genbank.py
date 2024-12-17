"""
Count the kmers in a genbank file and print them  out in order
"""


import os
import sys
import argparse
from itertools import product
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Count the kmers in a fasta file')
parser.add_argument('-g', help='genbank file', required=True)
parser.add_argument('-m', help='method (do all or just count existing)', choices=['all', 'count'], default='count')
parser.add_argument('-s', help='K-mer size, (default=11)', type=int, default=11)
parser.add_argument('-n', help='number of kmers to print (default = 10)', type=int, default=10)
parser.add_argument('-r', help='also count kmers on reverse complement', action='store_true')
parser.add_argument('-v', help='verbose output', action='store_true')
args = parser.parse_args()

gbk = SeqIO.parse(args.g, 'genbank')
count = {}

if args.m == 'all':
    for record in gbk:
        rcseq = None
        if args.r:
            rcseq = record.seq.reverse_complement()
        for k in product("ATGC", repeat=args.s):
            sk = "".join(k)
            count[sk] = record.seq.upper().count(sk)
            if rcseq:
                count[sk] += rcseq.upper().count(sk)
elif args.m == 'count':
    for record in gbk:
        posn = 0
        while posn < len(record.seq) - args.s:
            count[str(record.seq[posn:posn + args.s])] = count.get(record.seq[posn:posn + args.s], 0) + 1
            posn += 1
        if args.r:
            rcseq = record.seq.reverse_complement()
            posn = 0
            while posn < len(rcseq) - args.s:
                count[str(rcseq[posn:posn + args.s])] = count.get(rcseq[posn:posn + args.s], 0) + 1
                posn += 1



sorted_counts = sorted(count.items(), key=lambda item: item[1], reverse=True)
if args.n > len(sorted_counts):
    if args.v:
        print(f"Outputting all {len(sorted_counts)} kmers", file=sys.stderr)
    args.n = len(sorted_counts)
for i in range(args.n):
    print(sorted_counts[i])