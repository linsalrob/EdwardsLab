"""
A really simple fastq trimmer. You should probably use prinseq or trimmomatic
"""

import os
import sys
import argparse
from roblib import message, stream_fastq

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='fastq file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-w', '--forward', help="forward primer")
    parser.add_argument('-r', '--reverse', help='reverse primer (note we do not reverse complement, we look at the right end]')
    parser.add_argument('--maxfwd', type=int, default=1e6,
                        help='do not trim more than these bp from the start (does not include primer length)')
    parser.add_argument('--maxrev', type=int, default=1e6,
                        help='do not trim more than these bp from the end (does not include primer length)')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not args.forward and not args.reverse:
        message("Either --forward or --reverse primer must be specified otherwise nothing will be removed")
        sys.exit(-1)

    fwd = args.forward.upper()
    rev = args.reverse.upper()
    
    with open(args.o, 'w') as out:
        for sid, seqid, seq, qual in stream_fastq(args.f):
            original = [seq, qual]
            trimmed = False
            if fwd and fwd in seq.upper():
                idx = seq.upper().index(fwd)
                if idx < args.maxfwd:
                    if idx > 10:
                        message(f"WARNING: Trimming forward primer {fwd} from {sid} starting at position {idx}", "PINK")
                    seq = seq[idx+len(args.forward):]
                    qual = qual[idx+len(args.forward):]
                    trimmed = True
            if rev and rev in seq.upper():
                idx = seq.upper().rindex(rev)
                if idx < args.maxrev:
                    if idx > 10:
                        message(f"WARNING: Trimming reverse primer {rev} from {sid} starting at rightward position {idx}", "PINK")
                    seq = seq[0:idx]
                    qual = qual[0:idx]
                    trimmed = True
            if trimmed and args.v:
                message(f"Trimmed {sid}\n{original[0]}\n{seq}", "BLUE")
            out.write(f"@{seqid}\n{seq}\n+\n{qual}\n")
