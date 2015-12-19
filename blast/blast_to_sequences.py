import os
import sys


import argparse

import rob

__author__ = 'Rob Edwards'


def read_blast_file(filename, query=True, evalue=10, bitscore=0):
    """
    Read the blast output file and return a dict of hits that has contig, start, stop.
# crAssphage_C    NODE_1_length_14386_cov_54.5706_ID_1703 94.64   1157    62      0       82      1238    4975    3819    0.0     1794    1238    14386

    Using -outfmt '6 std qlen slen' the columns are:
        0: query
        1: database
        2: percent id
        3: alignment length
        4: gaps
        5: mismatches
        6: query start
        7: query end
        8: database start
        9; database end
        10: e value
        11: bit score
        12: query len
        13: subject len

    :param query: Retrieve hits from the query sequence (if false, we'll get them from the database sequence)
    :type query: bool
    :param bitscore: minimum bitscore to be included as a hit
    :type bitscore: int
    :param evalue: maximum E value to be included as a hit
    :type evalue: float
    :param filename: blast output filename (in tab separated text format)
    :type filename: str
    :return: dictionary of contigs, starts and stops for all hits
    :rtype: dict
    """

    if not os.path.exists(filename):
        sys.exit("{} not found\n".format(filename))



    hits = {}
    with open(filename, 'r') as fin:
        for l in fin:
            p = l.strip().split("\t")
            for i in range(3, len(p)):
                if i == 2 or i == 10 or i == 11:
                    p[i] = float(p[i])
                else:
                    p[i] = int(p[i])
            if p[11] < bitscore:
                continue
            if p[10] > evalue:
                continue
            if query:
                hitname = p[1]
                contig, start, end = p[0], p[6], p[7]
            else:
                hitname = p[0]
                contig, start, end = p[1], p[8], p[9]
            if contig not in hits:
                hits[contig] = []
            rc = False
            if start > end:
                start, end, rc = end, start, True
                start -= 1
            else:
                start -= 1
            hits[contig].append((start, end, rc, hitname))
    return hits


def extract_sequences(fastafile, hits, addhitname=False):
    """
    Extract the sequences from a fasta file

    :param fastafile: The fasta file to get the sequences from
    :type fastafile: str
    :param hits: The dict of hits using contig, start, end
    :type hits: dict
    :return: A dict of the sequences with contig_start_end as ID and sequence as value
    :rtype: dict
    """

    sequences = {}
    if not os.path.exists(fastafile):
        sys.exit("{} not found\n".format(fastafile))

    fa = rob.read_fasta(fastafile)

    for contig in hits:
        if contig not in fa:
            sys.stderr.write("WARNING: {} was not found in {}\n".format(contig, fastafile))

        for tple in hits[contig]:
          seq = fa[contig][tple[0]:tple[1]]
          if tple[2]:
              seq = rob.rc(seq)
          loc = "_".join(map(str, [contig, tple[0]+1, tple[1]]))
          if addhitname:
              loc += " [hit={}]".format(tple[3])
          sequences[loc] = seq
    return sequences


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='extract sequences based on blast hits')
    parser.add_argument('-b', help='blast output file', required=True)
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-d', help='use database sequences (default: query sequences', action="store_true")
    parser.add_argument('-e', help='Maximum evalue (default = 10)', type=float)
    parser.add_argument('-s', help='Minimum bit score (default = 0)', type=float)
    parser.add_argument('-i', help='Add database (query) ID to output', action='store_true', default=False)
    args = parser.parse_args()

    usequery = not args.d
    useeval = 10
    usebits = 0
    if args.e:
        useeval = args.e
    if args.s:
        usebits = args.s
    blasthits = read_blast_file(args.b, query=usequery, evalue=useeval, bitscore=usebits)
    blastseqs = extract_sequences(args.f, blasthits, args.i)
    for i in blastseqs:
        print(">{}\n{}".format(i, blastseqs[i]))