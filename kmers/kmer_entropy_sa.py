"""
Count the kmers in a file and report entropy and eveness.
This version should run anywhere
"""

import os
import sys
import argparse
import gzip

from math import log2

BLUE = '\033[94m'
GREEN = '\033[92m'
ENDC = '\033[0m'
RED = '\033[91m'

def rc(dna):
    """
    Reverse complement a DNA sequence

    :param dna: The DNA sequence
    :type dna: str
    :return: The reverse complement of the DNA sequence
    :rtype: str
    """
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

def stream_fasta(fastafile,):
    """
    Stream a fasta file, one read at a time. Saves memory!

    :param fastafile: The fasta file to stream
    :type fastafile: str
    :return:The ID, and a single read
    :rtype:str, str
    """

    if not os.path.exists(fastafile):
        sys.stderr.write(f"{RED}FATAL: {fastafile} does not exist\n{ENDC}")
        sys.exit(2)

    try:
        if fastafile.endswith('.gz'):
            f = gzip.open(fastafile, 'rt')
        else:
            f = open(fastafile, 'r')
    except IOError as e:
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write("Message: \n" + str(e.message) + "\n")
        sys.exit("Unable to open file " + fastafile)

    posn = 0
    while f:
        # first line should start with >
        idline = f.readline()
        if not idline:
            break
        if not idline.startswith('>'):
            sys.exit("Do not have a fasta file at: {}".format(idline))

        idline = idline.strip().replace('>', '', 1)
        posn = f.tell()
        line = f.readline()
        seq = ""
        while not line.startswith('>'):
            seq += line.strip()
            posn = f.tell()
            line = f.readline()
            if not line:
                break
        f.seek(posn)
        yield idline, seq


def count_kmers(faf, k, verbose=False):
    """
    Count the kmers
    :param faf: fasta file
    :param k: kmer size
    :param verbose: more output
    :return: a dict of kmers
    """

    if verbose:
        sys.stderr.write(f"{GREEN}Counting kmers (k={k}) in {faf}{ENDC}\n")

    kmers = {}

    for id, seq in stream_fasta(faf):
        rcseq = rc(seq)
        posn = 0
        while posn < len(seq) - k - 1:
            kmers[seq[posn:posn+k]] = kmers.get(seq[posn:posn+k], 0) + 1
            kmers[rcseq[posn:posn + k]] = kmers.get(rcseq[posn:posn + k], 0) + 1
            posn += 1

    if verbose:
        sys.stderr.write(f"{BLUE}\tDone counting kmers (k={k}) in {faf}{ENDC}\n")

    return kmers

def shannon(kmers, verbose=False):
    """
    Calculate the shannon entropy
    :param kmers: the kmer dictionary
    :param verbose: more output
    :return: the shannon entropy of the kmers
    """

    if verbose:
        sys.stderr.write(f"{GREEN}Calculating Shannon's Entropy{ENDC}\n")
    t = sum(kmers.values())
    H = 0
    for x in kmers:
        H += (kmers[x] / t) * (log2(kmers[x]/t))

    return -H

def evenness(kmers, H=None, verbose=False):
    """
    Calculate the evenness
    :param kmers: the kmer dictionary
    :param H: shannon entropy (optional). If provided, we won't recalculate
    :param verbose: more output
    :return: the evenness of the kmers
    """

    if not H:
        H = shannon(kmers, verbose)
    S = len(kmers.keys())
    return H/log2(S)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count the kmers in a file and report entropy and eveness')
    parser.add_argument('-f', help='fasta file to count the entropy/evenness', required=True)
    parser.add_argument('-k', help='kmer size', required=True, type=int)
    parser.add_argument('-t', help='print field titles in output', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    kmers = count_kmers(args.f, args.k, args.v)
    H = shannon(kmers, args.v)
    e = evenness(kmers, H, args.v)

    if args.t:
        print("File\tK-mer size\tShannon's Entropy\tEvenness")
    print(f"{args.f}\t{args.k}\t{H}\t{e}")