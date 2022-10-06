"""
Heroically try and figure out the destiny of every read in a pair of fastq files!

Not sure it will work, but a bunch of parsing
"""

import os
import sys
import argparse
from roblib.sequences import stream_fastq_index

__author__ = 'Rob Edwards'

def seqids_idx(fqfile: str, verbose:bool = False) -> dict[str, int]:
    """
    Read a fastq file and return a dict of seq ids and locations
    """

    seq = {}
    if verbose:
        print(f"Indexing reads in {fqfile}", file=sys.stderr)
    for seqid, header, seq, qualscores, seq_position in stream_fastq_index(fqfile):
        seq[seqid] = seq_position
    return seq


def read_tsv(r1: str = None, r2: str = None, column:int = 0, verbose: bool = False) -> tuple[set, set]:
    """
    Find reads that are in the superfocus m8 files and return sets of r1 and r2
    Either sf1 or sf2 (or both) can be None
    :param r1: Reads 1 output tsv file
    :param r2: Reads 2 output tsv file
    :param column: column with the read ID
    :param verbose: more output
    """

    rs1 = set()
    rs2 = set()
    if r1:
        if verbose:
            print(f"Parsing superfocus output {r1}", file=sys.stderr)
        with open(r1, 'r') as f:
            for l in f:
                p = l.strip().split("\t")
                rs1.add(p[column])
    if r2:
        if verbose:
            print(f"Parsing superfocus output {r2}", file=sys.stderr)
        with open(r2, 'r') as f:
            for l in f:
                p = l.strip().split("\t")
                rs2.add(p[column])


    return rs1, rs2

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-1', '--r1', help='R1 fastq file', required=True)
    parser.add_argument('-2', '--r2', help='R2 fastq file', required=True)
    parser.add_argument('-l', '--unknown_r1', help='File to write unknown R1 fastqs', required=True)
    parser.add_argument('-r', '--unknown_r2', help='File to write unknown R2 fastqs', required=True)
    group = parser.add_argument_group('tsv files', 'Use this for each tsv file, but specify all three options for each file')
    group.add_argument('-t1', '--tsv_r1', help='tsv file reads 1', action='append')
    group.add_argument('-t2', '--tsv_r2', help='tsv file reads 2', action='append')
    group.add_argument('-c', '--tsv_col', help='tsv column', action='append')
    parser.add_argument('-s', '--summary', help='summary file to write reason for matching. R1 and R2 will be appended')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if len(args.tsv_r1) != len(args.tsv_r2) or len(args.tsv_r1) != len(args.tsv_col):
        print("Please provide the same number of options for -t1, -t2, and -c", file=sys.stderr)
        sys.exit(-1)

    summaryR1 = None
    summaryR2 = None
    if args.summary:
        summaryR1 = open("{args.summary}.R1", 'w')
        summaryR2 = open("{args.summary}.R2", 'w')

    reads1 = seqids_idx(args.r1, args.verbose)
    reads2 = seqids_idx(args.r2, args.verbose)
    unkr1 = set(reads1.keys())
    unkr2 = set(reads2.keys())
    for i,j in enumerate(args.tsv_r1):
        tr1, tr2 = read_tsv(args.tsv_r1[i], args.tsv_r2[i], args.tsv_col[i])
        unkr1.difference_update(tr1)
        unkr2.difference_update(tr2)
        if summaryR1:
            for r in tr1:
                summaryR1.write(f"{r}\t{args.tsv_r1[i]}")
            for r in tr2:
                summaryR2.write(f"{r}\t{args.tsv_r2[i]}")


    if summaryR1:
        summaryR1.close()
        summaryR2.close()

    with open(args.unknown_r1, 'w') as fout, open(args.r1, 'r') as fin:
        for u in unkr1:
            fin.seek(reads1[u])
            for i in range(4):
                fout.write(fin.readline())
    with open(args.unknown_r2, 'w') as fout, open(args.r2, 'r') as fin:
        for u in unkr2:
            fin.seek(reads2[u])
            for i in range(4):
                fout.write(fin.readline())
