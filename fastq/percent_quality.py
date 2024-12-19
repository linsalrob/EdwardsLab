"""
Calculate the percent of a read less than a quality score
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
    parser.add_argument('-q', help='quality score', type=int, required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    print("SeqID\tLength\tAverage Qual")
    for sid, seqid, seq, qual in stream_fastq(args.f):
        q2n =  list(qual_to_numbers(qual))
        less = [x for x in q2n if x < args.q]
        if len(q2n) == 0:
            av = 0
        else:
            av = len(less)/len(q2n)
        print(f"{sid}\t{len(seq)}\t{av}")



