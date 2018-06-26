"""
Find the mates of some reads
"""

import os
import sys
import argparse
import re
from roblib import stream_fastq

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
            base = re.sub('\.2$', '', p[0])
            if base not in reads:
                reads[base] = [None, None]
            reads[base][1] = p[0]

    return reads

def odd_reads(reads):
    """
    Report reads that don't have a mate with a blast match
    :param reads: the dict of reads from the blast file
    :return: two sets of reads that are missing, one from read 1 and one from read 2
    """

    m1 = set()
    m2 = set()
    for r in reads:
        if reads[r][0] and not reads[r][1]:
            sys.stderr.write("No read 2: {}\n".format(r))
            m2.add("{}.2".format(r))

        if not reads[r][0] and reads[r][1]:
            sys.stderr.write("No read 1: {}\n".format(r))
            m1.add("{}.1".format(r))

    return m1, m2


def print_reads(miss1, miss2, fq1, fq2):
    """
    Print the missing reads from the two fastq files
    :param miss1: the set of reads missing from fq1
    :param miss2: the set of reads missing from fq2
    :param fq1: the first fastq file
    :param fq2: the second fastq file
    :return:
    """

    for sid, allid, seq, qual in stream_fastq(fq1):
        if sid in miss1:
            sys.stdout.write("@{}\n{}\n+\n{}\n".format(allid, seq, qual))

    for sid, allid, seq, qual in stream_fastq(fq2):
        if sid in miss2:
            sys.stdout.write("@{}\n{}\n+\n{}\n".format(allid, seq, qual))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', help='blast output files. Please provide exactly two files', nargs=2)
    parser.add_argument('-q', help='fastq files. Please provide exactly two files', nargs=2)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    sys.stderr.write("Using blast files: {}\nUsing fastq files: {}\n".format(args.b, args.q))

    if len(args.b) != 2:
        sys.stderr.write("Please provide two, and exactly two, blast output files\n")
        sys.exit(-1)

    reads = find_mates(*args.b)
    miss1, miss2 = odd_reads(reads)
    print_reads(miss1, miss2, *args.q)
