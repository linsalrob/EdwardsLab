"""
Test streaming pairs of sequences but not a unit (nose) test. Sorry
"""

import os
import sys
import argparse

from roblib import bcolors
from roblib import stream_paired_fastq



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-l', help='R1 file', required=True)
    parser.add_argument('-r', help='R2 file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    for seqid, h1, s1, q1, h2, s2, q2 in stream_paired_fastq(args.l, args.r):
        print(f"{h1} :: {h2}")