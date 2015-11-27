import os
import sys
import argparse
import pysam
__author__ = 'Rob Edwards'

"""
Count the reads associated with one or more contigs.
"""


parser = argparse.ArgumentParser(description='Count reads to one or more contigs in a bam file')
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculate the coverage of each read in a bam file')
    parser.add_argument('-b', help='bam file to parse')
    parser.add_argument('-s', help='sam file to parse')
    parser.add_argument('-c', help='partial contig name. If provided will be matched using *is in*')
    args = parser.parse_args()

    samfile = None
    if args.b:
        samfile = pysam.AlignmentFile(args.b, "rb")
    elif args.s:
        samfile = pysam.AlignmentFile(args.s, "r")
    else:
        sys.exit("Either -s or -b must be specified")

    count = {}
    for read in samfile.fetch():
        # print("{} -> {} : {}".format(read.query_name, read.reference_name, read.query_alignment_length))
        if args.c and args.c not in read.reference_name:
            continue
        count[read.reference_name] = count.get(read.reference_name, 0) + 1

    for c in count:
        print("{}\t{}".format(c, count[c]))

