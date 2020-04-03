"""
Compare coding and non coding regions of the fasta file
"""

import os
import sys
import argparse
from roblib import bcolors, stream_fasta

def coding_versus_noncoding(sample, contigs, orfs, outputfile):
    """
    Count the number of coding vs. noncoding bases
    """
    sys.stderr.write(f"{bcolors.GREEN}Counting coding vs. non coding bases{bcolors.ENDC}\n")

    # read the DNA sequences
    dnalen = 0
    for seqid, seq in stream_fasta(contigs):
        dnalen += len(seq)

    # read the protein sequences
    protlen = 0
    for seqid, seq in stream_fasta(orfs):
        protlen += len(seq)

    with open(outputfile, 'w') as out:
        out.write(f"{sample}\tCoding vs non coding\t")
        out.write(f"[coding bp, total bp, fraction coding]\t")
        out.write(f"{protlen}\t{dnalen}\t{protlen / dnalen}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-s', help='sample name used in output', required=True)
    parser.add_argument('-c', help='contigs file of DNA sequences', required=True)
    parser.add_argument('-f', help='fasta file of dna sequences of genes', required=True)
    parser.add_argument('-o', help='output file', required=True)
    args = parser.parse_args()

    coding_versus_noncoding(args.s, args.c, args.f, args.o)