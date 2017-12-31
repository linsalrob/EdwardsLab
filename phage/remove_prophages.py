"""
Read some directories of files and generate fasta files from genomes that have the
prophage regions excisesd
"""

import os
import sys
import argparse
import roblib


def next_location(s, start_posn, delimiters):
    """
    Find the next instance of the delimiter in the list of delimiters in the string s
    :param s: a string to search through
    :param start_posn: current position to start at. Note we start at that posn. 0==start of string
    :param delimiters: the list of things that may be used to break that string
    :return: a tuple of the next instance or the length of the string if none of the delimiters found and the
             and the next place to start looking
    """

    next_match = len(s)
    enddelim = len(s)
    for d in delimiters:
        try:
            thisposn = s.index(d, start_posn)
        except ValueError:
            # this delimiter is not in the string
            continue
        if thisposn < next_match:
            next_match = thisposn
            enddelim = thisposn + len(d)

    return next_match, enddelim, s[start_posn:next_match]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate fasta files without prophages",
                                     epilog='Note that we expect a seed directory in -s for every phispy directory' +
                                            ' (in -p) and that the seed directory contains a contigs file')
    parser.add_argument('-p', help='phispy directory', required=True)
    parser.add_argument('-s', help='seed directory', required=True)
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    # read the directories
    for phispydir in os.listdir(args.p):
        if not os.path.exists(os.path.join(args.s, phispydir)):
            if args.v:
                sys.stderr.write(
                    'A seed directory matching the phispy directory {} was not found. Skipped\n'.format(phispydir))
            continue

        # read the phages
        phageseqs = []
        for phagefile in [x for x in os.listdir(os.path.join(args.p, phispydir)) if x.endswith('.fasta')]:
            for pid, phagecontig in roblib.stream_fasta(os.path.join(args.p, phispydir, phagefile)):
                phageseqs.append(phagecontig)

        # sort the phages longest to smallest
        phageseqs = sorted(phageseqs, key=len, reverse=True)

        # now read the sequences and split out on phageseqs
        if not os.path.exists(os.path.join(args.s, phispydir, "contigs")):
            sys.stderr.write(
                "Error: no contigs file was found at {}. Skipped\n".format(os.path.join(args.s, phispydir, "contigs")))

        os.mkdir(os.path.join(args.o, phispydir))
        out = open(os.path.join(args.o, phispydir, "contigs"), 'w')

        for gid, genomecontig in roblib.stream_fasta(os.path.join(args.s, phispydir, "contigs")):
            contigcount = 0
            posn = 0
            while posn < len(genomecontig):
                lowest, posn, ss = next_location(genomecontig, posn, phageseqs)
                contigcount += 1
                if len(ss) > 0:
                    out.write(">{}.{}\n{}".format(gid, contigcount, ss))
