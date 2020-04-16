"""
Filter the uniref 50 file using precalculated taxnonomy strings.

These are stored in ranks.tsv and are made by the program

python3 ~/GitHubs/EdwardsLab/ncbi/all_taxid_by_rank.py -t superkingdom -v > ranks.tsv

The uniref fasta file is split up using

gunzip -c uniref50.fasta.gz | awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100000==0){file=sprintf("uniref50_files/file%d.fa",n_seq);} print >> file; n_seq++; next;} { print >>
file; }'

We can then run this python code on the cluster to process all the samples in parallel
"""

import os
import sys
import argparse
import re
from roblib import stream_fasta, colours

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

def read_ranks(rf, verbose=False):
    """
    Read the ranks tsv file. This is just [taxid, rank]
    :param rf: ranks files
    :param verbose: more output
    :return: dict of ranks.
    """

    if verbose:
        sys.stderr.write(f"{colours.GREEN}Reading{rf}{colours.ENDC}\n")

    rank = {}
    with open(rf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            rank[p[0]] = p[1]
    return rank

def split_by_rank(faf, ranks, outdir, verbose=False):
    """
    Split the fasta file
    :param faf: fasta file
    :param ranks: dict of taxid and rank
    :param outdir: output directory
    :param verbose: more output
    :return:
    """

    s = re.compile('TaxID=(\d+)')

    if args.v:
        sys.stderr.write(f"{colours.GREEN}Splitting {faf}{colours.ENDC}\n")
    fhs = {}
    for seqid, seq in stream_fasta(faf):
        m = s.search(seqid)
        rnk = "root"
        if m:
            tid = m.groups()[0]
            if tid in ranks:
                rnk = ranks[tid]
        else:
            sys.stderr.write(f"{colours.RED}ERROR: No taxonomy in {seqid}{colours.ENDC}\n")

        if rnk not in fhs:
            fhs[rnk] = open(os.path.join(outdir, rnk + ".fasta"), 'w')
        fhs[rnk].write(f">{seqid}\n{seq}\n")

    for fh in fhs:
        fhs[fh].close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-r', help='ranks file (probably ranks.tsv)', required=True)
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    ranks = read_ranks(args.r, args.v)
    split_by_rank(args.f, ranks, args.o, args.v)