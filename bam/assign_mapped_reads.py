"""
For all the reads in a fastq file, figure out where they are in a series of
bam files
"""

import os
import re
import sys
import argparse
from roblib import stream_fastq
import pysam
from suffix_trees import STree


__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', '--fastq', help='input fastq file', required=True)
    parser.add_argument('-b', '--bam', help='bamfile(s). Probably need more than one!', action='append')
    parser.add_argument('-s', '--sample', help='sample name (used as column header)')
    parser.add_argument('-o', '--output', help='output file to write with all reads information')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.sample:
        sampleid = args.sample
    else:
        # figure out what the common prefix of our bam files is
        st = STree.STree(args.bam)
        sampleid = st.lcs()


    # read the fastq file and find all the reads
    reads = set()
    for seqid, header, seq, qualscores in stream_fastq(args.fastq):
        if seqid.endswith('.1') or seqid.endswith(r'\1') or seqid.endswith('.2') or seqid.endswith(r'\2'):
            seqid = seqid[:-2]
        reads.add(seqid)


    dests = {}
    # read each bam file
    for bamfile in args.bam:
        bamreader = pysam.AlignmentFile(bamfile, "rb")
        wanted = set()
        for read in bamreader.fetch(until_eof=True):
            seqid = read.query_name
            if seqid.endswith('.1') or seqid.endswith(r'\1') or seqid.endswith('.2') or seqid.endswith(r'\2'):
                seqid = seqid[:-2]
            if seqid in reads:
                dests[seqid] = bamfile
            else:
                sys.stderr.write(f"ERROR. Found sequence {seqid} in the bamfile that is not in the fastq file\n")

    counts = {x:0 for x in args.bam}
    counts['unassigned'] = 0
    counts['ambiguous'] = 0
    out = open(args.output, 'w')
    print(f"read\t{sampleid}", file=out)
    for r in reads:
        if r in dests:
            if len(dests[r]) > 1:
                counts['ambiguous'] = counts.get('ambiguous') + 1
                print(f"{r}\tamibguous", file=out)
            else:
                counts[dests[r]] = counts.get(dests[r]) + 1
                print(f"{r}\t{dests[r]}", file=out)
        else:
            counts['unassigned'] = counts.get('unassigned') + 1
            print(f"{r}\tunassigned", file=out)
    out.close()

    for c in counts:
        print(f"{c}\t{counts[c]}")

