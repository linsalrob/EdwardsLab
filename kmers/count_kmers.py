"""
Count the 11mers in a sequence
"""

import os
import sys
import argparse
from itertools import product
from roblib import read_fasta
from statistics import median


fa = read_fasta("/home/redwards/Desktop/83333.1.contigs")

for id in fa:
    count = []
    for k in product("ATGC", repeat=11):
        sk = "".join(k)
        count.append(fa[id].count(sk))

        print("id: {} len(seq): {} sum: {}  n: {} average: {} median: {} max: {}".format(
            id, len(fa[id]), sum(count), len(count), (1.0 * sum(count) / len(count)), median(count), max(count)
        ))


