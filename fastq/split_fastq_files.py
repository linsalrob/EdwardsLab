"""
Split paired end fastq files into different files
"""

import os
import sys
import argparse
from roblib import stream_fastq, message

__author__ = 'Rob Edwards'




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split PE fastq files into multiple files')
    parser.add_argument('-l', help='fastq file 1', required=True)
    parser.add_argument('-r', help='fastq file 2', required=True)
    parser.add_argument('-n', help='number of reads per file. Default = ', type=int, default=10000)
    parser.add_argument('-o', required=True,
                        help='stub of output file _n_R1.fastq and _n_R2.fastq will be added for files 1 to n')
    parser.add_argument('-d', help='Output directory. Default = fastq_split', default='fastq_split')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    os.makedirs(args.d, exist_ok=True)
    reads = {}
    filecounter = 1
    counter = 0
    out = open(os.path.join(args.d, f"{args.o}_{filecounter}.R1.fastq"))
    for sid, seqid, seq, qual in stream_fastq(args.l):
        counter += 1
        if counter > args.n:
            counter = 0
            filecounter += 1
            out.close()
            out = open(os.path.join(args.d, f"{args.o}_{filecounter}.R1.fastq"))
        out.write(f"@{seqid}\n{seq}\n+\n{qual}\n")
        reads[sid] = filecounter
    out.close()

    filecounter = 1
    counter = 0
    out = open(os.path.join(args.d, f"{args.o}_{filecounter}.R2.fastq"))
    for sid, seqid, seq, qual in stream_fastq(args.r):
        counter += 1
        if counter > args.n:
            counter = 0
            filecounter += 1
            out.close()
            out = open(os.path.join(args.d, f"{args.o}_{filecounter}.R2.fastq"))
        if sid in reads and reads[sid] != filecounter:
            message("ERROR: Different file locations for {sid}. Left read is in {reads[sid]}. Right read is in {filecounter}. Please ensure your sequences are paired!", "RED")
        if sid not in reads:
            message("ERROR: Found an unpaire read, {sid}")
        out.write(f"@{seqid}\n{seq}\n+\n{qual}\n")
        reads[sid] = filecounter
    out.close()


