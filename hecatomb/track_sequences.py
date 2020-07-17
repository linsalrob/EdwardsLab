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
    dna = {}
    qual = {}
    header = {}

    # initially didn't plan to keep all these :)
    for seqid, header, seq, qualscores in stream_fastq(args.f):
        dna[seqid] = seq.upper()
        qual[seqid] = qualscores
        header[seqid] = header

    changed = set()
    for step in range(1, 10):
        if args.v:
            message(f"Working on step {step}", "GREEN")
        fqf = os.path.join(args.q, f"step_{step}", f"{args.n}.s{step}.out.fastq")
        if not os.path.exists(fqf):
            message(f"FQ File {fqf} not found", "RED")
            continue
        seqs = []
        with open(os.path.join(args.o, f"step_{step}.text"), 'w') as out, open(os.path.join(args.o, f"step_{step}.fq"), 'w') as fqout:
            seen = set()
            for seqid, header, seq, qualscores in stream_fastq(fqf):
                seen.add(seqid)
                if seqid not in dna:
                    message(f"{seqid} is a different sequence id", "PINK")
                    continue
                if dna[seqid] != seq.upper():
                    out.write(f"{seqid}\n")
                    fqout.write(f"@{header[seqid]}\n{dna[seqid]}\n+\n{qual[seqid]}\n")
                    dna[seqid] = seq.upper()
                    changed.add(seqid)
            for seqid in dna:
                if seqid not in seen:
                    out.write(f"{seqid}\n")
                    fqout.write(f"@{header[seqid]}\n{dna[seqid]}\n+\n{qual[seqid]}\n")
                    changed.add(seqid)

    with open(os.path.join(args.o, "unchanged.txt"), 'w') as out:
        for s in dna:
            if s not in changed:
                out.write(f"{s}\n")





