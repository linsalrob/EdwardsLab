"""
Filter some reads based on a blast file.

Read the blast file, parse the list of queries based on e value or bit score and then
read a fastq file and filter those reads NOT in the blast file
"""

import os
import sys
import argparse
from roblib import stream_fastq, stream_blast_results

def filter_reads(blastf, eval, bitscore, verbose=False):
    """
    Filter our blast results
    :param blastf: the blast file
    :return: a set of the IDs to be filtered
    """

    blast = set()
    if verbose:
        sys.stderr.write(f"Parsing {blastf}\n")
    for br in stream_blast_results(blastf, verbose):
        if br.evalue > eval:
            continue
        if br.bitscore < bitscore:
            continue
        if verbose:
            sys.stderr.write("Found {}\n".format(br.query))
        blast.add(br.query)
    return blast

def read_fastq(fqfile, blast, verbose=False):
    """
    Read the fastq file and print only sequences we need
    :param fqfile:  The fastq file
    :param blast: the blast reads that matched (ie. reads to delete)
    :param verbose: more output
    :return:
    """

    for seqid, fullid, seq, qual in stream_fastq(fqfile):
        if seqid.startswith('@'):
            seqid = seqid[1:]
        if seqid in blast or fullid in blast:
            continue
        print("@{}\n{}\n+\n{}".format(fullid, seq, qual))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter reads from a blast file")
    parser.add_argument('-f', help='blast file', required=True)
    parser.add_argument('-q', help='fastq file', required=True)
    parser.add_argument('-e', help='maximum e value allowed (default=all)', type=float, default=1000)
    parser.add_argument('-b', help='minimum bit score allowed (default=all)', type=int, default=0)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    reads = filter_reads(args.f, args.e, args.b, args.v)
    read_fastq(args.q, reads)
