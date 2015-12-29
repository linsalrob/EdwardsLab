import os
import sys
import argparse
from roblib import read_fasta

__author__ = 'Rob Edwards'

parser = argparse.ArgumentParser(description='Convert a fasta file to fastq, faking the qual scores')
parser.add_argument('-f', help='fasta file', required=True)
parser.add_argument('-q', help='fastq output file', required=True)
parser.add_argument('-s', help='quality score. Default = 40', default=40, type=int)
args = parser.parse_args()

c = chr(args.s)

fa = read_fasta(args.f)
with open(args.q, 'w') as out:
    for i in fa:
        l = len(fa[i])
        out.write("@{}\n{}\n+\n{}\n".format(i, fa[i], l * c))

