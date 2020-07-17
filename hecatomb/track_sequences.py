"""
Track sequences through the hecatomb steps and identify those sequences that have changed
at each step.

This will allow us to make some test sets
"""

import os
import sys
import argparse
from roblib import stream_fastq, message

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Track hecatomb sequences")
    parser.add_argument('-f', help='initial fastq file', required=True)
    parser.add_argument('-n', help='base file name. Everything upto the _R1', required=True)
    parser.add_argument('-q', help='QC dir (default: %(default)s', default='QC')
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    os.makedirs(args.o, exist_ok=True)
    initial = {}
    for seqid, header, seq, qualscores in stream_fastq(args.f):
        initial[seqid] = seq.upper()

    for step in range(1, 10):
        fqf = os.path.join(args.q, f"step_{step}", f"{args.n}.s{step}.out.fastq")
        if not os.path.exists(fqf):
            message(f"FQ File {fqf} not found", "RED")
            continue
        seqs = []
        with open(os.path.join(args.o, f"step_{step}.text"), 'w') as out:
            for seqid, header, seq, qualscores in stream_fastq(fqf):
                if seqid not in initial:
                    message(f"{seqid} is a different sequence id", "PINK")
                    continue
                if initial[seqid] != seq.upper():
                    out.write(f"{seqid}\n")
                    initial[seqid] = seq.upper()







