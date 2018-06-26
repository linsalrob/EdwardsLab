"""
Find the mates of some reads
"""

import os
import sys
import argparse
import re

def find_mates(blastnf1, blastnf2):
    """
    Find our mates
    :param blastnf1: first blastn file
    :param blastnf2: second blastn file
    :return: a dict of reads mapping basename -> [read1, read2]
    """

    reads = {}
    with open(blastnf1, 'r') as fin:
        for l in fin:
            p = l.strip().split("\t")
            base = re.sub('\.1$', '', p[0])
            reads[base] = [p[0], None]

    with open(blastnf2, 'r') as fin:
        for l in fin:
            p = l.strip().split("\t")
            base = re.sub('\.1$', '', p[0])
            if base not in reads:
                reads[base] = [None, None]
            reads[base][1] = p[0]

    return reads

def odd_reads(reads):
    """
    Report reads that don't have a mate with a blast match
    :param reads: the dict of reads from the blast file
    :return:
    """

    for r in reads:
        if reads[1] and not reads[2]:
            sys.stderr.write("No read 2: {}\n".format(r))

        if not reads[1] and reads[2]:
            sys.stderr.write("No read 1: {}\n".format(r))





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', help='blast output files (use two -b)', action='append')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    if len(args.b) != 2:
        sys.stderr.write("Please provide two, and exactly two, blast output files\n")
        sys.exit(-1)

    reads = find_mates(*args.b)