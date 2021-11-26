"""
Given a fastq file with barcode=barcode01 in the defline, split the fastqs into separate files.

Wrote to be standalone so don't need to clone all of github
"""

import os
import re
import sys
import argparse
import gzip


__author__ = 'Rob Edwards'

class FastqFormatError(Exception):
    """
    Exception raised for sequences not being paired properly.

    :param message: explanation of the error
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


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

    os.makedirs(args.d, exist_ok=True)
    srch = re.compile(r'barcode=(\S+)')
    files = {}

    for seqid, header, seq, qualscores in stream_fastq(args.f):
        mt = srch.findall(header)
        if not mt:
            sys.stderr.write(f"Error parsing {header}\n")
            continue
        bc = mt[0]
        if bc not in files:
            files[bc] = open(os.path.join(args.d, f"{bc}.fastq"), "w")
        files[bc].write(f"@{header}\n{seq}\n+\n{qualscores}\n")
    for f in files:
        files[f].close()
