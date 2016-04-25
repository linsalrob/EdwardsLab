"""
Given a BAM file and a set of fastq file(s) pull all the reads
(including pairs) that are in the BAM file. Note that if only
one of the pairs is in the BAM file we pull the other pair.

Assumes that the reads end .1 and .2 for paired end
"""

import os
import sys

import argparse
import pysam
import re
import rob

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam to fastq')
    parser.add_argument('-b', help='bam file', required=True)
    parser.add_argument('-f', help='fastq file(s) ... one or more can be specified', required=True, action='append')
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.b, "rb")
    wanted = set()
    for read in bamfile.fetch(until_eof=True):
        wanted.add(read.query_name)
        t = read.query_name
        if t.endswith('.1'):
            t = re.sub('\.1$', '.2', t)
            wanted.add(t)
        elif t.endswith('.2'):
            t = re.sub('\.2$', '.1', t)
            wanted.add(t)
    for f in args.f:
        for fq in rob.stream_fastq(f):
            if fq[0] in wanted:
                print("@" + fq[1] + "\n" + fq[2] + "\n+\n" + fq[3])

