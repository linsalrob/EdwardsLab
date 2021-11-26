"""
Given a fastq file with barcode=barcode01 in the defline, split the fastqs into separate files.

Wrote to be standalone so don't need to clone all of github
"""

import os
import re
import sys
import argparse

__author__ = 'Rob Edwards'


import os
import sys
import gzip

import subprocess
from .rob_error import SequencePairError, FastqFormatError
from .colours import colours, message

__author__ = 'Rob Edwards'


def read_fasta(fname: str, whole_id: bool = True, qual: bool = False) -> dict:
    """
    Read a fasta file and return a hash.

    If wholeId is set to false only the first part of the ID
    (upto the first white space) is returned

    :param fname: The file name to read
    :param whole_id: Whether to keep the whole id, or trim to first whitespace (default = all)
    :param qual: these are quality scores (so add a space between lines!)
    :return: dict
    """

    try:
        if fname.endswith('.gz'):
            f = gzip.open(fname, 'rt')
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
            if qual:
                seq += " " + line
            else:
                seq += line

    seqs[seqid] = seq.strip()
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
            raise FastqFormatError(f"The file does not appear to be a four-line fastq file at line {linecounter}")
        header = header.strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seqid = seqid.replace('@', '')
        seq = qin.readline().strip()
        linecounter += 1
        qualheader = qin.readline()
        if not qualheader.startswith("+"):
            raise FastqFormatError(f"The file does not appear to be a four-line fastq file at line {linecounter}")
        linecounter += 1
        qualscores = qin.readline().strip()
        linecounter += 1
        header = header.replace('@', '', 1)
        if len(qualscores) != len(seq):
            raise FastqFormatError(f"The sequence and qual scores are not the same length at line {linecounter}")
        yield seqid, header, seq, qualscores



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input fastq file', required=True)
    parser.add_argument('-d', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    srch = re.compile(r'barcode=(\S+)')
    files = {}
    for seqid, header, seq, qualscores in stream_fastq(args.f):
        mt = srch.findall(header)
        if not mt:
            sys.stderr.write(f"Error parsing {header}\n")
            continue
        bc = mt[0]
        if bc not in files:
            files[bc] = open(os.path.join(args.d, f"{bc}.fastq"))
        bc.write(f"@{header}\n{seq}\n+\n{qualscores}\n")
    for f in files:
        files[f].close()
