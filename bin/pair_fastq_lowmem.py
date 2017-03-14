"""
An alternate method to pair fastq files
"""

import os
import sys

import argparse

import gzip


def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores.

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """

    if fqfile.endswith('.gz'):
        qin = gzip.open(fqfile, 'rb')
    else:
        qin = open(fqfile, 'r')

    while True:
        header = qin.readline()
        if not header:
            break
        header = header.strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seq = qin.readline()
        seq = seq.strip()
        qualheader = qin.readline()
        qualscores = qin.readline()
        qualscores = qualscores.strip()
        header = header.replace('@', '', 1)
        yield seqid, header, seq, qualscores


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Pair fastq files, writing all the pairs to separate files and the unmapped reads to separate files")
    parser.add_argument('-l', help='Pair #1 reads file')
    parser.add_argument('-r', help='Pair #2 reads file')
    args = parser.parse_args()

    # read the first file into a data structure
    seqs = []
    for (seqid, header, seq, qual) in stream_fastq(args.l):
        seqid = seqid.replace('.1', '')
        seqs.append(seqid)

    rp = open("{}.paired.fastq".format(args.r.replace('.fastq', '').replace('.gz', '')), 'w')
    ru = open("{}.singles.fastq".format(args.r.replace('.fastq', '').replace('.gz', '')), 'w')

    # read the first file into a data structure
    seen = set()
    for (seqid, header, seq, qual) in stream_fastq(args.r):
        seqid = seqid.replace('.2', '')
        if seqid in seqs:
            seqs.remove(seqid)
            rp.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")
        else:
            ru.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")

    ru.close()
    rp.close()

    lp = open("{}.paired.fastq".format(args.l.replace('.fastq', '').replace('.gz', '')), 'w')
    lu = open("{}.singles.fastq".format(args.l.replace('.fastq', '').replace('.gz', '')), 'w')

    for (seqid, header, seq, qual) in stream_fastq(args.l):
        seqid = seqid.replace('.1', '')
        if seqid not in seqs:
            # we have seen this so we deleted it
            lp.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")
        else:
            lu.write("@" + header + "\n" + seq + "\n+\n" + qual + "\n")

    lp.close()
    lu.close()