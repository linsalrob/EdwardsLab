"""
Pair up the two fastq files using a bloom filter instead of a dict. This should be a memory efficient
solution.
"""

import os
import sys

import argparse
from bitstring import BitArray
import gzip


def line_to_id(line):
    """
    Convert a line to its ID
    :param line: The line from the fastq file
    :type line: str
    :return: The sequence ID without the .1 or .2
    :rtype: str
    """
    line = line.strip()
    seqid = line.split(' ')[0]
    # here we reverse the string, do the first replace of 1. (reverse!), and rereverse the string
    if seqid.endswith('.1'):
        seqid = seqid[::-1].replace('1.', '', 1)[::-1]
    elif seqid.endswith('_1'):
        seqid = seqid[::-1].replace('1_', '', 1)[::-1]
    elif seqid.endswith('.2'):
        seqid = seqid[::-1].replace('2.', '', 1)[::-1]
    elif seqid.endswith('_2'):
        seqid = seqid[::-1].replace('2_', '', 1)[::-1]
    else:
        sys.exit("We can not figure out whether " + l + " is a .1 or a .2 read")

    return seqid


def id_to_index(seqid):
    """
    Convert a string to a poisitive integer
    :param seqid: The sequence ID
    :type seqid: str
    :return: The positive integer corresponding to the ID
    :rtype: int
    """

    return hash(seqid) & 0x7FFFFFFF



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Pair fastq files, writing all the pairs to separate files and the unmapped reads to separate files")
    parser.add_argument('-l', help='Pair #1 reads file', required=True)
    parser.add_argument('-r', help='Pair #2 reads file', required=True)
    parser.add_argument('-n', help='Number of sequences in the file. This is used to determine the size of the bloom filter. The default is 500,000,000', type=int, default=500000000)
    args = parser.parse_args()

    # define the bloom filter to be of size N. This should be at least 1.3x the number of reads.
    sys.stderr.write("Initializing the bloom filter\n")
    lbf = BitArray(args.n)
    lbf.set(0)
    sys.stderr.write("Completed initializing the bloom filter\n")

    counter = 0

    # read and populate the bloom filter
    if args.l.endswith('.gz'):
        qin = gzip.open(args.l, 'rb')
    else:
        qin = open(args.l, 'r')

    while True:
        l = qin.readline()
        if not l:
            break
        counter += 1
        seqid = line_to_id(l)
        idx = id_to_index(seqid) % args.n
        lbf.set(1, idx)
        seq = qin.readline()
        plus = qin.readline()
        qual = qin.readline()

    # warn if we have too many sequences
    if counter * 1.3 > args.n:
        sys.stderr.write("Warning, we have {} sequences, but we only used a bloom filter of {}. There may be some false positives. Try adjusting -n\n".format(counter, args.n))


    # set up a bloom filter for the unpaired reads. We assume that only ~10% of the reads will be unpaired.
    rbf = BitArray(args.n / 10)
    rbf.set(0)

    rp = open("{}.paired.fastq".format(args.r.replace('.fastq', '').replace('.gz', '')), 'w')
    ru = open("{}.singles.fastq".format(args.r.replace('.fastq', '').replace('.gz', '')), 'w')

    # now read and write the right file
    if args.r.endswith('.gz'):
        qin = gzip.open(args.r, 'rb')
    else:
        qin = open(args.r, 'r')
    while True:
        l = qin.readline()
        if not l:
            break
        seq = qin.readline()
        plus = qin.readline()
        qual = qin.readline()
        seqid = line_to_id(l)
        idx = id_to_index(seqid) % args.n
        if lbf[idx]:
            rp.write("{}{}{}{}".format(l, seq, plus, qual))
        else:
            ru.write("{}{}{}{}".format(l, seq, plus, qual))

        rbf.set(1, idx)
    rp.close()
    ru.close()

    # finally, read and write the left file

    lp = open("{}.paired.fastq".format(args.l.replace('.fastq', '').replace('.gz', '')), 'w')
    lu = open("{}.singles.fastq".format(args.l.replace('.fastq', '').replace('.gz', '')), 'w')

    if args.l.endswith('.gz'):
        qin = gzip.open(args.l, 'rb')
    else:
        qin = open(args.l, 'r')
    while True:
        l = qin.readline()
        if not l:
            break
        seq = qin.readline()
        plus = qin.readline()
        qual = qin.readline()
        seqid = line_to_id(l)
        idx = id_to_index(seqid) % args.n
        if rbf[idx]:
            lp.write("{}{}{}{}".format(l, seq, plus, qual))
        else:
            lu.write("{}{}{}{}".format(l, seq, plus, qual))

    lp.close()
    lu.close()