import os
import sys
import gzip

import subprocess

__author__ = 'Rob Edwards'


def read_fasta(fname, whole_id=True):
    """
    Read a fasta file and return a hash.

    If wholeId is set to false only the first part of the ID
    (upto the first white space) is returned
    """

    try:
        if fname.endswith('.gz'):
            f = gzip.open(fname, 'rb')
        elif fname.endswith('.lrz'):
            f = subprocess.Popen(['/usr/bin/lrunzip', '-q', '-d', '-f', '-o-', fname], stdout=subprocess.PIPE).stdout
        else:
            f = open(fname, 'r')
    except IOError as e:
        sys.stderr.write(str(e) + "\n")

        sys.stderr.write("Message: \n" + str(e.message) + "\n")

        sys.exit("Unable to open file " + fname)

    seqs = {}
    seq = ""
    seqid = ""
    for line in f:
        line = line.rstrip('\r\n')
        if line.startswith(">"):
            if seqid != "":
                seqs[seqid] = seq
                seq = ""
            seqid = line.replace(">", "", 1)
            if not whole_id and seqid.count(" ") > 0:
                seqids = seqid.split(" ")
                seqid = seqids[0]
        else:
            seq += line

    seqs[seqid] = seq
    return seqs


def readFasta(file, whole_id=True):
    """
    Read a fasta file and return a hash.

    If wholeId is set to false only the first part of the ID (upto the first white space) is returned
    """
    return read_fasta(file, whole_id)


def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores.

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """

    with open(fqfile, 'r') as qin:
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


