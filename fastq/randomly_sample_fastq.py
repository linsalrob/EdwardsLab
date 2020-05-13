"""
Randomly sample a fraction of a fastq file and just create one file
"""

import os
import sys
import argparse
from random import sample
from roblib import stream_fastq

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Randomly sample a single fastq file")
    parser.add_argument('-f', help='fastq file to sample', required=True)
    parser.add_argument('-o', help='output file name', required=True)
    parser.add_argument('-p', help='percent of the file to sample', required=True, type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    sequences = []
    for seqid, header, seq, qualscores in stream_fastq(args.f):
        sequences.append([header, seq, qualscores])

    with open(args.o, 'w') as out:
        for s in sample(sequences, args.p):
            out.write(f"@{s[0]}\n{s[1]}\n+\n{s[2]}\n")
