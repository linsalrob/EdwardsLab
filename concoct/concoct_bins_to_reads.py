"""
Extract the reads associated with bins in a concoct output

We need to make several connections:
    1. bins to contigs
    2. contigs to reads
    3. find both mates of the reads
    4. write those to a fastq file for the bin

"""

import os
import sys
import argparse

from roblib import bcolors
import pysam

from roblib import bcolors, stream_paired_fastq, stream_fastq


def contig_clusters(cct, verbose=False):
    """
    Read the concoct output file and return a hash of contigs->clusters
    :param cct:
    :param verbose:
    :return:
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Reading {cct}{bcolors.ENDC}\n")
    clusters = {}
    with open(cct, 'r') as f:
        for l in f:
            p = l.strip().split(",")
            clusters[p[0]] = p[1]
    return clusters

def contig_reads(bamdir, clusters, verbose=False):
    """
    Read the directory of bam files and return a list of read ids that map to clusters
    :param bamdir: the directory of bam files
    :param clusters: hash of [contigid, cluster]
    :param verbose: more output
    :return: a dict of [readid, cluster]
    """


    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Converting bams to read clusters{bcolors.ENDC}\n")

    reads = {}
    for bamf in os.listdir(bamdir):
        if not bamf.endswith('.bam'):
            continue

        bam = pysam.AlignmentFile(os.path.join(bamdir, bamf), 'rb')
        for c in clusters:
            for rds in bam.fetch(c):
                if rds.qname in reads:
                    sys.stderr.write(f"{bcolors.PINK}WARNING: {rds.qname} mapped to two different locations{bcolors.ENDC}\n")
                    continue
                reads[rds.qname] = clusters[c]
    return reads

def write_sequences(reads, outdir, leftfq, rightfq, singlefq = None, verbose = False):
    """
    Write the sequences out to a file
    :param reads: the dict of reads and bins
    :param outdir: the output dir to write to
    :param leftfq: the left reads
    :param rightfq: the right reads
    :param singlefq: the single reads (optional)
    :param verbose: more output
    :return:
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Writing sequences\n{bcolors.ENDC}")

    if not os.path.exists(outdir):
        os.mkdir(outdir)


    files = {}
    psc = 0
    for seqid, header1, seq1, qualscores1, header2, seq2, qualscores2 in stream_paired_fastq(leftfq, rightfq):
        if seqid in reads:
            if reads[seqid] not in files:
                files[reads[seqid]] = [
                    open(os.path.join(outdir, reads[seqid] + ".R1.fastq"), 'w'),
                    open(os.path.join(outdir, reads[seqid] + ".R2.fastq"), 'w')
                ]
            files[reads[seqid]][0].write(f"@{header1}\n{seq1}\n+\n{qualscores1}\n")
            files[reads[seqid]][1].write(f"@{header2}\n{seq2}\n+\n{qualscores2}\n")
            psc += 1

    singlefiles = {}
    sc = 0
    if singlefq:
        for seqid, header, seq, qualscores in stream_fastq(singlefq):
            if seqid in reads:
                if reads[seqid] not in singlefiles:
                    singlefiles[reads[seqid]] = open(os.path.join(outdir, reads[seqid] + ".single.fastq"), 'w')
            singlefiles[reads[seqid]].write(f"@{header}\n{seq}\n+\n{qualscores}\n")
            sc += 1

    for f in files:
        files[f][0].close()
        files[f][1].close()
    for f in singlefiles:
        singlefiles[f].close()

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Wrote {psc} paired end sequences and {sc} single reads\n{bcolors.ENDC}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Concoct to reads")
    parser.add_argument('-c', help='concoct merged clusters file', required=True)
    parser.add_argument('-b', help='directory of indexed bam files', required=True)
    parser.add_argument('-l', help='left reads', required=True)
    parser.add_argument('-r', help='right reads', required=True)
    parser.add_argument('-s', help='singleton reads (optional)')
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    clusters = contig_clusters(args.c, args.v)
    rds = contig_reads(args.b, clusters, args.v)
    write_sequences(rds, args.o, args.l, args.r, args.s, args.v)
