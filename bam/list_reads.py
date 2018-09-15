"""
List all the reads that map to a bam file
"""

import os
import sys
import argparse
import pysam

parser = argparse.ArgumentParser(description="List all the reads in a bam file")
parser.add_argument('-b', help='bam file', required=True)
parser.add_argument('-v', help='verbose output', action="store_true")
args = parser.parse_args()

bamfile = pysam.AlignmentFile(args.b, "rb")
for read in bamfile.fetch(until_eof=True):
    print(read.query_name)
