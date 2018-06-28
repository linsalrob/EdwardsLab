import os
import sys
import gzip

import subprocess

__author__ = 'Rob Edwards'




def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores.

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """

    needsdecode = False
    if fqfile.endswith('.gz'):
        qin = gzip.open(fqfile, 'rb')
        needsdecode = True
    else:
        qin = open(fqfile, 'r')

    while True:
        header = qin.readline()
        if needsdecode:
            header = header.decode('UTF-8')
        if not header:
            break
        header = header.strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seq = qin.readline()
        if needsdecode:
            seq = seq.decode('UTF-8')
        seq = seq.strip()
        qualheader = qin.readline()
        qualscores = qin.readline()
        if needsdecode:
            qualscores = qualscores.decode('UTF-8')
        qualscores = qualscores.strip()
        header = header.replace('@', '', 1)
        yield seqid, header, seq, qualscores


