"""
Convert embl format to a couple of other formats
"""

import os
import sys

import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert EMBL format sequence files")
    parser.add_argument('-e', help='EMBL input file', required=True)
    parser.add_argument('-f', help='Store fasta file', action="store_true")
    parser.add_argument('-g', help='Store genbank file', action="store_true")
    parser.add_argument('-o', help='output file base name')
    args = parser.parse_args()

    outfile = args.e
    if args.o:
        outfile = args.o

    sin =  SeqIO.read(args.e, 'embl')

    if args.f:
        outfafile = outfile + ".fasta"
        sout = SeqIO.write(sin, outfafile, 'fasta')
        print("Wrote {} records to a fasta file".format(sout))

    if args.g:
        outgbfile = outfile + ".gbk"
        sout = SeqIO.write(sin, outgbfile, 'genbank')
        print("Wrote {} records to a genbank file".format(sout))
