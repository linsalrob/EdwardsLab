"""
Compare the same fq files in two different directories
"""

import os
import sys
import argparse

from roblib import stream_fastq, colours

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count the fastq sequences for the same files in two directories and compare')
    parser.add_argument('-d', '--dir1', help='first directory of fastq files', required=True)
    parser.add_argument('-e', '--dir2', help='second directory of fastq files', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    print(f"Filename\t# seqs in {args.dir1}\t# bp in {args.dir1}\t# seqs in {args.dir2}\tbp in {args.dir2}")
    for fq in set(os.listdir(args.dir1)).union(set(os.listdir(args.dir2))):
        if args.v:
            sys.stderr.write(f"{colours.GREEN}Reading: {fq}{colours.ENDC}")
        data = [fq,0,0,0,0]
        if os.path.exists(os.path.join(args.dir1, fq)):
            for seqid, header, seq, qualscores in stream_fastq(os.path.join(args.dir1, fq)):
                data[1] += 1
                data[2] += len(seq)
        if os.path.exists(os.path.join(args.dir2, fq)):
            for seqid, header, seq, qualscores in stream_fastq(os.path.join(args.dir2, fq)):
                data[3] += 1
                data[4] += len(seq)
        print("\t".join(map(str, data)))