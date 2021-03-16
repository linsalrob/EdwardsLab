import os
import sys
import argparse
from roblib import read_fasta
from random import randint, shuffle

__author__ = 'Rob Edwards'

parser = argparse.ArgumentParser(description='Convert a fasta file to fastq, faking the qual scores. Either -f or -n is required')
parser.add_argument('-f', help='fasta file')
parser.add_argument('-n', help='number of sequences', type=int)
parser.add_argument('-q', help='fastq output file', required=True)
parser.add_argument('-s', help='quality score. Default = 40', default=40, type=int)
parser.add_argument('-l', help='length of sequences for artificial sequences. Default = 150', default=150, type=int)
parser.add_argument('-r', help='random quality scores between 5 and 40', action='store_true')
args = parser.parse_args()

if not args.f and not args.n:
    sys.stderr.write("FATAL: Either -f or -n is required. Please provide one or the other")
    sys.exit(0)

c = chr(args.s)

if args.f:
    fa = read_fasta(args.f)
    with open(args.q, 'w') as out:
        for i in fa:
            l = len(fa[i])
            q = l * c
            if args.r:
                q=""
                for s in range(l):
                    q = q + chr(randint(33, 125))
            out.write("@{}\n{}\n+\n{}\n".format(i, fa[i], q))
    exit(0)

if args.n:
    bases = ['A', 'T', 'G', 'C']
    with open(args.q, 'w') as out:
        for i in range(args.n):
            seq = ''
            qual = ''
            for j in range(args.l):
                shuffle(bases)
                seq = seq + bases[0]
                qual = qual + chr(randint(33, 125))
            out.write("@FakeSeq{}\n{}\n+\n{}\n".format(i, seq, qual))
