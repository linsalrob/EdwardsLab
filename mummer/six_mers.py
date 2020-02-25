"""
Create a list of all six mers and filter for those that are also  in the list of reverse complement 6 mers
"""

import os
import sys
import argparse

from roblib import bcolors
from roblib import rc

from itertools import product

def countit():
    kmer = 6
    bases = {'A', 'C', 'T', 'G'}

    fwd = set()
    rev = set()
    for s in product(bases, repeat = kmer):
        seq = "".join(s)
        fwd.add(seq)
        rev.add(rc(seq))

    count = 0
    for f in fwd:
        if f not in rev:
            print(f)
        else:
            count+=1

    print(f"Checked {count}")

if __name__ == '__main__':
    countit()