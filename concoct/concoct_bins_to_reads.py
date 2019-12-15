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
    allbamfs = os.listdir(bamdir)
    totalf = len(allbamfs)
    fc = 0
    if verbose:
        sys.stderr.write(f"{bcolors.BLUE}")
    for bamf in os.listdir(bamdir):
        if not bamf.endswith('.bam'):
            continue
        fc += 1
        if verbose:
            sys.stderr.write(f"\tProcessing file {fc} of {totalf}\r")
            sys.stderr.flush()
        bam = pysam.AlignmentFile(os.path.join(bamdir, bamf), 'rb')
        # find the reference names
        for s in bam.get_index_statistics():
            if s.contig not in clusters:
                continue
            for rds in bam.fetch(s.contig):
                if rds.qname not in reads:
                    reads[rds.qname] = {clusters[s.contig]}
                else:
                    reads[rds.qname].add(clusters[s.contig])

    if verbose:
        # finish the writing
        sys.stderr.write(f"{bcolors.ENDC}\n")

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
            for clst in reads[seqid]:
                if clst not in files:
                    files[clst] = [
                        open(os.path.join(outdir, clst + ".R1.fastq"), 'w'),
                        open(os.path.join(outdir, clst + ".R2.fastq"), 'w')
                    ]
                files[clst][0].write(f"@{header1}\n{seq1}\n+\n{qualscores1}\n")
                files[clst][1].write(f"@{header2}\n{seq2}\n+\n{qualscores2}\n")
            psc += 1

    singlefiles = {}
    sc = 0
    if singlefq:
        for seqid, header, seq, qualscores in stream_fastq(singlefq):
            if seqid in reads:
                for clst in reads[seqid]:
                    if clst not in singlefiles:
                        singlefiles[clst] = open(os.path.join(outdir, clst + ".single.fastq"), 'w')
                singlefiles[clst].write(f"@{header}\n{seq}\n+\n{qualscores}\n")
            sc += 1

    for f in files:
        files[f][0].close()
        files[f][1].close()
    for f in singlefiles:
        singlefiles[f].close()

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Wrote {psc} paired end sequences and {sc} single reads\n{bcolors.ENDC}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert concoct to reads",
        epilog="The concoct clusters file is ususally called clustering_merged.csv.")
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
