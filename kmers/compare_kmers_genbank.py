"""
given a set of files, calculate the number of k-mers shared between each file
"""

import os
import sys
import argparse
from itertools import combinations
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Count the kmers in a fasta file')
parser.add_argument('-g', help='genbank files (multiple required)', required=True, action='append')
parser.add_argument('-s', help='K-mer size, (default=11)', type=int, default=11)
args = parser.parse_args()

count = {}
for gbkfile in args.g:
    gbk = SeqIO.parse(gbkfile, 'genbank')
    count[gbkfile] = set()
    for record in gbk:
        posn = 0
        while posn < len(record.seq) - args.s:
            count[gbkfile].add(str(record.seq[posn:posn + args.s]))
            posn += 1

for tple in combinations(count.keys(), 2):
    common = count[tple[0]].intersection(count[tple[1]])
    union = count[tple[0]].union(count[tple[1]])
    print(f"{tple[0]}\t{tple[1]}\t{len(common)/len(union)}")
