import argparse
import os
import sys
from rob import stream_fastq

__author__ = 'Rob Edwards'

"""
Make sure that the pairs are really pairs, and write the singletons to a file
"""



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check paired end files and make sure the pairs match up')
    parser.add_argument('-l', help='The file where the reads end /1 (the option is a lowercase L)', required=True)
    parser.add_argument('-r', help='The file where the reads end /2', required=True)
    args = parser.parse_args()

    lseq = {}
    for (seqid, header, seq, qual) in stream_fastq(args.l):
        if not seqid.endswith('/1'):
            sys.stderr.write("Sequence {} in {} does not appear to be a read /1\n".format(seqid, args.l))
            continue
        seqid = seqid.replace('/1', '')
        lseq[seqid] = [header, seq, qual]

    rseq = {}
    for (seqid, header, seq, qual) in stream_fastq(args.r):
        if not seqid.endswith('/2'):
            sys.stderr.write("Sequence {} in {} does not appear to be a read /2\n".format(seqid, args.l))
            continue
    seqid = seqid.replace('/2', '')
    rseq[seqid] = [header, seq, qual]

    singles = set()
    with open("{}.paired".format(args.l), 'w') as lout:
        with open("{}.paired".format(args.r), 'w') as rout:
            for seqid in lseq:
                if seqid in rseq:
                    lout.write("@{}\n{}\n+\n{}\n".format(*lseq[seqid]))
                    rout.write("@{}\n{}\n+\n{}\n".format(*rseq[seqid]))
                else:
                    singles.add(seqid)

    with open("{}.singles".format(args.l), 'w') as lout:
        for seqid in singles:
            lout.write("@{}\n{}\n+\n{}\n".format(*lseq[seqid]))

    with open("{}.singles".format(args.r), 'w') as rout:
        for seqid in rseq:
            if seqid not in lseq:
                rout.write("@{}\n{}\n+\n{}\n".format(*rseq[seqid]))