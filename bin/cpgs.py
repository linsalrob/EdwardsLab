"""
Read a fastq file and count CpGs
"""

import os
import sys
import argparse

from roblib import bcolors
from roblib import stream_fastq

def countcpgs(fqfile):
    """
    Count the CpGs in a file
    :param fqfile: the fastq file
    :return:
    """

    count = {}
    for seqid, header, seq, qual in stream_fastq(fqfile):
        cg = seq.count('CG')
        count[cg] = count.get(cg, 0) + 1
    return count


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count CGs in a fastq file')
    parser.add_argument('-f', help='fastq file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    count = countcpgs(args.f)
    for c in sorted(list(count.keys())):
        print(f"{c}\t{count[c]}")
