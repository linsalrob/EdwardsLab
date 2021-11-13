"""
Given the output of correlation_clustering (a json file), and a fasta file, we write the bins to a directory

Optional: we can take an index file that maps contigs in the correlation file to IDs in the fasta file
"""
import json
import os
import sys
import argparse
import gzip

__author__ = 'Rob Edwards'



def read_fasta(fname: str, verbose: bool) -> dict:
    """
    Read a fasta file and return a dict of sequences

    :param fname: The file name to read
    :param verbose: more output
    :return: dict
    """
    if verbose:
        sys.stderr.write(f"Reading fasta file {fname}\n")

    try:
        if fname.endswith('.gz'):
            f = gzip.open(fname, 'rt')
        else:
            f = open(fname, 'r')
    except IOError as e:
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write("Message: \n" + str(e.message) + "\n")
        sys.exit("Unable to open file " + fname)

    seqs = {}
    seq = ""
    seqid = ""
    for line in f:
        line = line.rstrip('\r\n')
        if line.startswith(">"):
            if seqid != "":
                seqs[seqid] = seq.strip()
                seq = ""
            seqid = line.replace(">", "", 1)
        else:
            seq += line

    seqs[seqid] = seq.strip()
    return seqs

def read_idx(idxf: str, verbose: bool) -> dict:
    """
    Read the index file and return a dict of numbers -> sequence IDs
    :param idxf: index file
    :param verbose: more output
    :return: dict of nos->ids
    """
    if verbose:
        sys.stderr.write(f"Reading index file {idxf}\n")

    idx = {}
    with open(idxf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            idx[p[0]] = p[1]
    return idx

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', '--fasta', help='fasta file of DNA sequences (contigs)', required=True)
    parser.add_argument('-c', '--clusters', required=True,
                        help='cluster json file. The output from correlation_clustering.py')
    parser.add_argument('-d', '--directory', help='Output directory', required=True)
    parser.add_argument('-i', '--idx', help='Optional index file that maps numbers in clusters to fasta sequence IDs')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    idx2id = None
    if args.idx:
        idx2id = read_idx(args.idx, args.verbose)

    seqs = read_fasta(args.fasta, args.verbose)
    with open(args.clusters, 'r') as binf:
        bins = json.load(binf)

    os.makedirs(args.directory, exist_ok=True)

    for bin in bins['clusters']:
        with open(os.path.join(args.directory, f"atavide_bin_{bin}.fa"), 'w') as out:
            for sid in bins['clusters'][bin]:
                if idx2id:
                    try:
                        seqid = idx2id[sid]
                    except KeyError as e:
                        sys.stderr.write(f"ERROR: {sid} was not found in your sequence ID file {args.idx}\n")
                        sys.exit(2)
                else:
                    seqid = sid
                if seqid not in seqs:
                    sys.stderr.write(f"ERROR: {seqid} was not found in your fasta file {args.fasta}\n")
                    sys.exit(2)
                out.write(f">{seqid}\n{seqs[seqid]}\n")


