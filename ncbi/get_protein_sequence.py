"""
Get a set of protein sequences from NCBI
"""

import os
import sys

import argparse

from Bio import Entrez, SeqIO
from time import sleep
from random import randint

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get a set of protein sequences from NCBI")
    parser.add_argument('-f', help='File of IDs to get', required=True)
    parser.add_argument('-o', help='Output file', required=True)
    args = parser.parse_args()

    # retrieve a GI number from GenBank
    Entrez.email = 'raedwards@gmail.com'  # set this so NCBI knows who to complain to

    out = open(args.o, 'w')
    with open(args.f, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            handle = Entrez.efetch(db="protein", id=p[0], rettype="gbwithparts", retmode="text")
            out.write(handle.read())
            sleep(randint(0, 5))
    out.close()
