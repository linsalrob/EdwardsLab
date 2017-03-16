"""
Convert a fastq file to a fasta file. Note in this case I just ignore the quailty scores!
"""

import os
import sys
import argparse
from roblib import stream_fastq

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input fastq file', required=True)
    parser.add_argument('-o', help='output fasta file', required=True)
    args = parser.parse_args()

    with open(args.o, 'w') as out:
        for (sid, label, seq, qual) in stream_fastq(args.f):
            out.write(">{}\n{}\n".format(sid, seq))
