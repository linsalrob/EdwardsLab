"""
Count the kmers in a file and report entropy and eveness

NOTE: DO NOT USE THIS VERSION

It is designed to work with an old version of python (3.3) and you should not use it!
"""

import os
import sys
import argparse

from itertools import product

import json
from math import log2

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




def stream_fasta(fastafile, whole_id=True):
    """
    Stream a fasta file, one read at a time. Saves memory!

    :param fastafile: The fasta file to stream
    :type fastafile: str
    :param whole_id: Whether to return the whole id (default) or just up to the first white space
    :type whole_id:bool
    :return:A single read
    :rtype:str, str
    """

    try:
        if fastafile.endswith('.gz'):
            f = gzip.open(fastafile, 'rt')
        elif fastafile.endswith('.lrz'):
            f = subprocess.Popen(['/usr/bin/lrunzip', '-q', '-d', '-f', '-o-', fastafile], stdout=subprocess.PIPE).stdout
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
        if not whole_id:
            idline = idline.split(" ")[0]
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



def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores.

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """

    if fqfile.endswith('.gz'):
        qin = gzip.open(fqfile, 'rt')
    else:
        qin = open(fqfile, 'r')

    linecounter = 0
    while True:
        header = qin.readline()
        linecounter += 1
        if not header:
            break
        if not header.startswith("@"):
            raise FastqFormatError("The file does not appear to be a four-line fastq file at line ")
        header = header.strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seqid = seqid.replace('@', '')
        seq = qin.readline().strip()
        linecounter += 1
        qualheader = qin.readline()
        if not qualheader.startswith("+"):
            raise FastqFormatError("The file does not appear to be a four-line fastq file at ")
        linecounter += 1
        qualscores = qin.readline().strip()
        linecounter += 1
        header = header.replace('@', '', 1)
        if len(qualscores) != len(seq):
            raise FastqFormatError("The sequence and qual scores are not the same length at line")
        yield seqid, header, seq, qualscores



def count_kmers(faf, type, k, jsonout=None, verbose=False):
    """
    Count the kmers
    :param faf: fasta file
    :param type: str either fasta or fastq
    :param k: kmer size
    :param verbose: more output
    :return: a dict of kmers
    """

    if verbose:
        sys.stderr.write("Counting kmers (k={}) in {}\n".format(k, faf))

    kmers = {}

    if type == "fasta":
        for id, seq in stream_fasta(faf):
            rcseq = rc(seq)
            posn = 0
            while posn < len(seq) - k - 1:
                kmers[seq[posn:posn+k]] = kmers.get(seq[posn:posn+k], 0) + 1
                kmers[rcseq[posn:posn + k]] = kmers.get(rcseq[posn:posn + k], 0) + 1
                posn += 1

    if type == "fastq":
        for id, fullid, seq, qual in stream_fastq(faf):
            rcseq = rc(seq)
            posn = 0
            while posn < len(seq) - k - 1:
                kmers[seq[posn:posn+k]] = kmers.get(seq[posn:posn+k], 0) + 1
                kmers[rcseq[posn:posn + k]] = kmers.get(rcseq[posn:posn + k], 0) + 1
                posn += 1

    if jsonout:
        with open(jsonout, 'w') as out:
            json.dump({faf : kmers}, out)

    if verbose:
        sys.stderr.write("\tDone counting kmers\n")

    return kmers

def shannon(kmers, verbose=False):
    """
    Calculate the shannon entropy
    :param kmers: the kmer dictionary
    :param verbose: more output
    :return: the shannon entropy of the kmers
    """

    if verbose:
        sys.stderr.write("Calculating Shannon's Entropy\n")
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
    parser.add_argument('-f', help='fasta file to count the entropy/evenness')
    parser.add_argument('-q', help='fastq file to count the entropy/evenness')
    parser.add_argument('-k', help='kmer size', required=True, type=int)
    parser.add_argument('-j', help='json output for kmer counts')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.f:
        kmers = count_kmers(args.f, 'fasta', args.k, args.j, args.v)
    elif args.q:
        kmers = count_kmers(args.q, 'fastq', args.k, args.j, args.v)
    else:
        sys.stderr.write("FATAL: Please supply either a fasta file or a fastq file\n")
        sys.exit(-1)

    H = shannon(kmers, args.v)
    e = evenness(kmers, H, args.v)

    if args.f:
        print("{}\t{}\t{}\t{}".format(args.f, args.k, H, e))
    else:
        print("{}\t{}\t{}\t{}".format(args.q, args.k, H, e))
