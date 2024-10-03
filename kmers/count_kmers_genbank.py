"""
Count the kmers in a genbank file
"""


import os
import sys
import argparse
from itertools import product
from roblib import median
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Count the kmers in a fasta file')
parser.add_argument('-g', help='genbank file', required=True)
parser.add_argument('-s', help='K-mer size, (default=11)', type=int, default=11)
args = parser.parse_args()

gbk = SeqIO.parse(args.g, 'genbank')
for record in gbk:
    count = []
    for k in product("ATGC", repeat=args.s):
        sk = "".join(k)
        count.append(record.seq.upper().count(sk))

    print("id: {} len(seq): {} sum: {}  n: {} average: {} median: {} max: {}".format(
        id, len(record.seq), sum(count), len(count), (1.0 * sum(count) / len(count)), median(count), max(count)
    ))


