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
            base = p[0][:p[0].rindex(".1")]
            reads[base] = [p[0], None]

    with open(blastnf2, 'r') as fin:
        for l in fin:
            p = l.strip().split("\t")
            base = p[0][:p[0].rindex(".2")]
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

    miss = set()
    for r in reads:
        if reads[r][0] and not reads[r][1]:
            sys.stderr.write("No read 2: {}\n".format(r))
            miss.add(r)

        if not reads[r][0] and reads[r][1]:
            sys.stderr.write("No read 1: {}\n".format(r))
            miss.add(r)

    return miss


def print_reads(miss, fq1, fq2):
    """
    Print the missing reads from the two fastq files
    :param miss1: the set of reads missing from fq1
    :param miss2: the set of reads missing from fq2
    :param fq1: the first fastq file
    :param fq2: the second fastq file
    :return:
    """


    bn = re.search('/(\w+)_pass_1.fastq', fq1)
    if not bn:
        sys.stderr.write(f"Can't parse the base filename from {fq1}\n")
        sys.exit(-1)

    fqo1 = bn.groups()[0] + "_missed_1.fastq"
    fqo2 = bn.groups()[0] + "_missed_2.fastq"
    if os.path.exists(fqo1):
        sys.stderr.write(f"Not overwrting {fqo1}\n")
        sys.exit(-1)

    if os.path.exists(fqo2):
        sys.stderr.write(f"Not overwrting {fqo2}\n")
        sys.exit(-1)

    with open(fqo1, 'w') as out:
        sys.stderr.write("Finding reads from {}\n".format(fq1))
        c =  0
        for sid, allid, seq, qual in stream_fastq(fq1):
            c += 1
            if not c % 100000:
                sys.stderr.write(".")
                sys.stderr.flush()
            test = sid[:sid.rindex(".1")].replace('@', '', 1)
            if test in miss:
                out.write("@{}\n{}\n+\n{}\n".format(allid, seq, qual))
                out.flush()

    with open(fqo2, 'w') as out:
        sys.stderr.write("\nFinding reads from {}\n".format(fq2))
        c=0
        for sid, allid, seq, qual in stream_fastq(fq2):
            c += 1
            if not c % 100000:
                sys.stderr.write(".")
                sys.stderr.flush()

            test = sid[:sid.rindex(".2")].replace('@', '', 1)
            if test in miss:
                out.write("@{}\n{}\n+\n{}\n".format(allid, seq, qual))
                out.flush()
    sys.stderr.write("\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', help='blast output files. Please provide exactly two files', nargs=2)
    parser.add_argument('-q', help='fastq files. Please provide exactly two files. If these are not supplied we print the unique reads and exit', nargs=2)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    sys.stderr.write("Using blast files: {}\nUsing fastq files: {}\n".format(args.b, args.q))

    if len(args.b) != 2:
        sys.stderr.write("Please provide two, and exactly two, blast output files\n")
        sys.exit(-1)

    reads = find_mates(*args.b)
    miss = odd_reads(reads)
    if args.q:
        print_reads(miss, *args.q)
