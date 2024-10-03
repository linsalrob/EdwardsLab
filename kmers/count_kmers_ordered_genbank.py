"""
Count the kmers in a genbank file and print them  out in order
"""


import os
import sys
import argparse
from itertools import product
from roblib import median
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Count the kmers in a fasta file')
parser.add_argument('-g', help='genbank file', required=True)
parser.add_argument('-m', help='method (do all or just count existing)', choices=['all', 'count'], default='count')
parser.add_argument('-s', help='K-mer size, (default=11)', type=int, default=11)
parser.add_argument('-n', help='number of kmers to print (default = 10)', type=int, default=10)
args = parser.parse_args()

gbk = SeqIO.parse(args.g, 'genbank')
count = {}
if args.m == 'all':
    for record in gbk:
        for k in product("ATGC", repeat=args.s):
            sk = "".join(k)
            count[sk] = record.seq.upper().count(sk)
elif args.m == 'count':
    for record in gbk:
        posn = 0
        while posn < len(record.seq) - args.s:
            count[str(record.seq[posn:posn + args.s])] = count.get(record.seq[posn:posn + args.s], 0) + 1
            posn += 1

sorted_counts = sorted(count.items(), key=lambda item: item[1], reverse=True)
for i in range(args.n):
    print(sorted_counts[i])