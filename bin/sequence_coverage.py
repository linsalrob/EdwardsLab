import os
import sys

__author__ = 'Rob Edwards'

"""
Read a BAM file and calculate the coverage for each read
"""

import pysam
import argparse






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculate the coverage of each read in a bam file')
    parser.add_argument('-b', help='bam file to parse')
    parser.add_argument('-s', help='sam file to parse')
    args = parser.parse_args()

    samfile = None
    if args.b:
        samfile = pysam.AlignmentFile(args.b, "rb")
    elif args.s:
        samfile = pysam.AlignmentFile(args.s, "r")
    else:
        sys.exit("Either -s or -b must be specified")


    for read in samfile.fetch():
        print("{}".format(read))