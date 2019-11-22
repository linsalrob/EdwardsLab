"""
Create all possible pairwise permutations of a fasta file so we can run
it through needleman wunsch. There are probably better ways to do this
like incorporating the nw code into this code, but this is easy. 
It will just make a really big fasta file!
"""

import os
import sys
import argparse
from roblib import bcolors
from roblib import read_fasta
from itertools import combinations
__author__ = 'Rob Edwards'

def write_permutations(faf, outputf, verbose=False):
    """
    Create and write all the permutations
    """

    fa = read_fasta(faf, whole_id=False)
    ids = list(fa.keys())

    with open(outputf, 'w') as out:
        for tple in combinations(ids, 2):
            out.write(f">{tple[0]}\n{fa[tple[0]]}\n>{tple[1]}\n{fa[tple[1]]}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    write_permutations(args.f, args.o, args.v)



