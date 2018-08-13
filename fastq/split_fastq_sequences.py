"""
Break fastq sequences into smaller chunks. Eg. from nanopore reads we want smaller pieces
"""

import os
import sys
import argparse
from roblib import stream_fastq

__author__ = 'Rob Edwards'

def rewrite_fastq(inf, outf, sz, verbose):
    """
    Rewrite a fastq file
    :param inf: input fastq file
    :param outf: output fastq file
    :param sz: size of the DNA sequences to write
    :param verbose: more output
    :return:
    """

    with open(outf, 'w') as out:
        seqcounter = 0
        for seqid, header, seq, qual in stream_fastq(inf):
            posn = 0
            while (posn < len(seq)-sz):
                seqcounter += 1
                out.write("@{} {}\n{}\n+\n{}\n".format(seqcounter, header, seq[posn:posn+sz], qual[posn:posn+sz]))
                posn += sz
            if posn  < len(seq):
                out.write("@{} {}\n{}\n+\n{}\n".format(seqcounter, header, seq[posn:], qual[posn:]))





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read a fastq file and write the sequences smaller. Note does not preserve IDs!')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-s', help='DNA fragment size', required=True, type=int)
    parser.add_argument('-o', help='output file name', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    rewrite_fastq(args.f, args.o, args.s, args.v)