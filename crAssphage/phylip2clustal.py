"""
A simple BioPython converter to move from phylip to clustal formats for the alignments
"""

import os
import sys

import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert an alignment file from phylip format to clustal format")
    parser.add_argument('-i', help='Alignment input file', required=True)
    parser.add_argument('-o', help='Output file name (optional: default = input file with clustal appended)')
    args = parser.parse_args()

    outfile = args.i + ".clustal"
    if args.o:
        outfile = args.o

    records=SeqIO.parse(args.i, 'phylip')

    with open(outfile, 'w') as out:
        SeqIO.write(records, out, 'clustal')