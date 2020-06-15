"""
Filter based on quality/length
"""

import os
import sys
import argparse

from roblib import stream_fastq, qual_to_numbers

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='file', required=True)
    parser.add_argument('-q', help='minimum average quality', type=float, default=0)
    parser.add_argument('-l', help='minimum length', type=int, default=0)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    print("SeqID\tLength\tAverage Qual")
    for sid, seqid, seq, qual in stream_fastq(args.f):
        if len(seq) <= args.l:
            continue
        q2n =  list(qual_to_numbers(qual))
        av = sum(q2n)/len(q2n)
        if av > args.q:
            print(f"{sid}\t{len(seq)}\t{av}")