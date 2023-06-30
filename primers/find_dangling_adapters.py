"""
Find adapter sequences at the end of reads
"""

import os
import sys
import argparse
from roblib import stream_fastq

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input fastq file', required=True)
    parser.add_argument('-a', help='adapter sequence. Default = AGATCGGAAGAGCACACGTCTGAACTCCAGTCA')
    parser.add_argument('-m', help="minimum length of sequence", type=int, default=5)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.a:
        adapter = args.a
    else:
        adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

    for sid, allid, seq, qual in stream_fastq(args.f):
        end = len(adapter)
        while (end >= args.m):
            if adapter[0:end] in seq:
                offset = seq.find(adapter[0:end]) + end
                print(f"{sid}\t{adapter[0:end]}\t{offset}\t{len(seq) - offset}")
                break
            end -= 1
