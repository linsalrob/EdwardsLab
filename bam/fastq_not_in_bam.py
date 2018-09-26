"""
Extract the reads from a fastq file that are not in the bam file
"""

import os
import sys
import argparse
import pysam
from roblib import stream_fastq

def reads_from_bam(bamf, verbose=False):
    """
    Read the reads in a bam file
    :param bamf: bam file
    :param verbose: more output
    :return: a set of read ids in the file
    """

    reads = set()
    bamfile = pysam.AlignmentFile(bamf, "rb")
    for read in bamfile.fetch(until_eof=True):
        reads.add(read.query_name)
    if verbose:
        sys.stderr.write("There are {} reads in the bamfile\n".format(len(reads)))
    return reads

def extract_fastq(fqf, reads, verbose):
    """
    Extract the reads from the fastq file
    :param fqf: fastq file
    :param reads: set of reads to ignore
    :param verbose: more output
    :return:  nada
    """

    for (sid, label, seq, qual) in stream_fastq(fqf):
        if sid.startswith('@'):
            sid = sid[1:]
        if sid not in reads:
            if verbose:
                sys.stderr.write("Keeping: {}  -->  {}\n".format(sid, label))
            print("@{}\n{}\n+\n{}".format(label, seq, qual))
        elif verbose:
            sys.stderr.write("Skipping: {}  -->  {}\n".format(sid, label))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract fastq reads that are not in the bamfile")
    parser.add_argument('-f', help='fastq file', required=True)
    parser.add_argument('-b', help='bamfile', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    reads = reads_from_bam(args.b, args.v)
    extract_fastq(args.f, reads, args.v)