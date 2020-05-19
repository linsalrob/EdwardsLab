"""
Convert a fasta/quality files to a fastq file. I can't believe I'm writing this in 2020
"""

import os
import sys
import argparse

from roblib import read_fasta, write_fastq, message

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-q', help='quality file', required=True)
    parser.add_argument('-o', help='output fastq file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.f) and not os.path.exists(args.q):
        message("FATAL: either {args.f} or {args.q} not found", "RED")
        sys.exit(-1)

    fa = read_fasta(args.f, True, False) 
    qu = read_fasta(args.q, True, True)

    write_fastq(fa, qu, args.o, args.v)
