"""
Model the frequency that we find different kmers in DNA sequences
"""

import os
import sys
import argparse
from random import randint
from itertools import product
from statistics import median

__author__ = 'Rob Edwards'

def generate_random_seq(length):
    """
    Generate a random sequence of length len
    :param length: the length to generate
    :return: str
    """

    bases = {1: "A", 2: "G", 3: "C", 4: "T"}
    seq=""
    for i in range(length):
        seq += bases[randint(1,4)]
    return seq



#seq = generate_random_seq(1000000)
seq = generate_random_seq(4194304)
count = []
for k in product("ATGC", repeat=11):
    sk="".join(k)
    count.append(seq.count(sk))



print("sum: {}  n: {} average: {} median: {}".format(
    sum(count), len(count), (1.0*sum(count)/len(count)), median(count)
))

