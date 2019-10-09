import argparse
import os
import sys

import pysam
from roblib import read_fasta, bcolors

__author__ = 'Rob Edwards'

def list2ids(listf, verbose=False):
    """ Read a list of ids, one per line
    :param listf: file of ids
    :param verbose: more output
    :return: a set of IDS
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN} Reading IDs from text file:  {listf}{bcolors.ENDC}")
    i=set()
    with open(listf, 'r') as f:
        for l in f:
            i.add(l.strip().split(" ")[0])
    return i

def fasta2ids(faf, verbose=False):
    """
    Extract IDs from a fasta file
    :param faf: fasta file
    :param verbose: more output
    :return: a set of IDS
    """
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN} Reading IDs from fasta file: {faf}{bcolors.ENDC}")

    f = read_fasta(faf, whole_id=False)
    return set(f.keys())

def qual2fastq(quals, verbose=False):
    """
    Convert a list of quality scores to a single fastq line

    :param quals: A list of quality scores
    :type quals: list
    :return: A fastq quality string
    :rtype: str
    """
    quality = [chr(q + 33) for q in quals]
    return "".join(quality)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam to fastq')
    parser.add_argument('-b', help='bam file', required=True)
    parser.add_argument('-l', help='file to output fastq left reads', required=True)
    parser.add_argument('-r', help='file to output fastq right reads', required=True)
    parser.add_argument('-i', help='list of ids (one per line) to keep')
    parser.add_argument('-f', help='fasta file of sequences whose ids to keep')
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.b, "rb")
    i = None
    if args.l:
        i = list2ids(args.l, args.v)
    if args.f:
        if i:
            i.update(fasta2ids(args.f, args.v))
        else:
            i = fasta2ids(args.f, args.v)

    if args.v:
        sys.stderr.write(f"{bcolors.BLUE} Writing to {args.l} and {args.r}{bcolors.ENDC}")

    with open(args.l, 'w') as left:
        with open(args.r, 'w') as right:
            for read in bamfile.fetch(until_eof=True):
                if i and read not in i:
                    continue
                # construct the output string
                if read.query_qualities:
                    ostr = "@{}\n{}\n+\n{}".format(read.query_name, read.query_sequence, qual2fastq(read.query_qualities))
                else:
                    ostr = "@{}\n{}\n+\n".format(read.query_name, read.query_sequence)


                if read.endswith("1") or read.endswith("l"):
                    left.write(ostr)
                else:
                    right.write(ostr)


