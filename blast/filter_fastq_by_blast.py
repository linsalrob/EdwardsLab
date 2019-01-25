#!/usr/bin/env python3

"""
Filter a fastq file from a blast results file. The fastq file must be the query sequence.
"""

import os
import sys
import argparse
from roblib import stream_blast_results, stream_fastq, bcolors

__author__ = 'Rob Edwards'

def read_blast(blastf, maxeval, minlen, minid, verbose=False):
    """
    Read the blast file and only keep those matches that exceed our parameters.
    :param blastf: blast file to read
    :param maxeval: maximum ID to keep
    :param minlen: minimum length to keep
    :param minid: minimum % id to keep
    :param verbose: more output
    :return: a set of matches
    """

    results = set()

    for b in stream_blast_results(blastf, verbose=verbose):
        if b.alignment_length < minlen:
            continue
        if b.evalue > maxeval:
            continue
        if b.percent_id < minid:
            continue

        results.add(b.query)
    return results

def filter_fastq(fqf, br, matchout=None, nomatchout=None, verbose=False):
    """
    Filter the fastq file and print out matches or no matches
    :param fqf: The fastq file to filter
    :param br: the set of query blast results
    :param matchout: The file to write matches to
    :param nomatchout: the file to write no matches to
    :param verbose: more output
    :return: nothing
    """

    mo = open(matchout, 'w')
    nmo = open(nomatchout, 'w')

    matches = 0
    nonmatches = 0
    for sid, allid, seq, qual in stream_fastq(fqf):
        if sid in br:
            if matchout:
                mo.write(f"@{allid}\n{seq}\n+\n{qual}\n")
            matches += 1
        else:
            if nomatchout:
                nmo.write(f"@{allid}\n{seq}\n+\n{qual}\n")
                nonmatches += 1
    sys.stderr.write(f"{bcolors.GREEN}FINISHED:{bcolors.ENDC} Sequences Matched: {matches} Sequences without match {nonmatches}\n")


if __name__ == "__main__":
    maxeval = 1e-10
    minlen  = 50
    minid   = 75
    parser = argparse.ArgumentParser(description='Filter fastq files based on blast results')
    parser.add_argument('-f', help='fastq file to filter', required=True)
    parser.add_argument('-b', help='blast output file (using -outfmt 6 std)', required=True)
    parser.add_argument('-m', help='file to write the sequences that match the blast file to')
    parser.add_argument('-n', help='file to write the sequences that DO NOT match the blast file to')
    parser.add_argument('-e', help='Maximum E value cut off (default={})'.format(maxeval), default=maxeval)
    parser.add_argument('-l', help='Minimum alignment length cut off (default={})'.format(minlen), default=minlen)
    parser.add_argument('-i', help='Minimum percent id cut off (default={})'.format(minid), default=minid)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not args.m and not args.n:
        sys.stderr.write(f"{bcolors.FAIL}ERROR: Either -m or -n must be specified{bcolors.ENDC}\n")
        sys.exit(-1)

    b = read_blast(args.b, args.e, args.l, args.i, args.v)
    if not args.m:
        args.m = None
    if not args.n:
        args.n = None
    filter_fastq(args.f, b, args.m, args.n, args.v)
