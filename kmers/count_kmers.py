"""
Count the 11mers in a sequence
"""

import os
import sys
import argparse
from itertools import product
from roblib import read_fasta
from statistics import median

parser = argparse.ArgumentParser(description='Count the kmers in a fasta file')
parser.add_argument('-f', help='fasta file', required=True)
parser.add_argument('-s', help='K-mer size, (default=11)', type=int, size=11)
args = parser.parse_args()

fa = read_fasta(args.f)

for id in fa:
    count = []
    for k in product("ATGC", repeat=args.s):
        sk = "".join(k)
        count.append(fa[id].count(sk))

        print("id: {} len(seq): {} sum: {}  n: {} average: {} median: {} max: {}".format(
            id, len(fa[id]), sum(count), len(count), (1.0 * sum(count) / len(count)), median(count), max(count)
        ))


