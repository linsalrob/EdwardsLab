"""
Use the prophage tbl file to acurrately and correctly write the prophages from a set of predictions
"""

import os
import sys
import argparse
from operator import itemgetter
from roblib import stream_fasta
import re


def parse_prophage_tbl(phispydir):
    """
    Parse the prophage table and return a dict of objects

    :param phispydir: The phispy directory to find the results
    :return: dict
    """
    if not os.path.exists(os.path.join(phispydir, "prophage.tbl")):
        sys.stderr.write("FATAL: The file prophage.tbl does not exist\n")
        sys.stderr.write("Please run create_prophage_tbl.py -d {}\n".format(phispydir))
        sys.exit(-1)

    p = re.compile('^(.*)_(\d+)_(\d+)$')

    locations = {}
    with open(os.path.join(phispydir, "prophage.tbl"), 'r') as f:
        for l in f:
            (ppid, location) = l.strip().split("\t")
            m = p.search(location)
            (contig, beg, end) = m.groups()
            beg = int(beg)
            end = int(end)
            if beg > end:
                (beg, end) = (end, beg)
            if contig not in locations:
                locations[contig] = []
            locations[contig].append((beg, end, ppid))
    return locations


def parse_contigs(locations, gdir, odir):
    """
    Parse the contigs file and print non-prophage regions

    :param locations: the locations hash from the phispy directory
    :param gdir: the genome directory that contains the contigs file
    :param odir: the output directory
    :return: None
    """

    out = open(os.path.join(odir, "pp.fasta"), 'w')
    for contig, seq in stream_fasta(os.path.join(gdir, "contigs")):
        if contig not in locations:
            continue
        ses = sorted(locations[contig], key=itemgetter(0))
        for start, end, ppid in ses:
            out.write(">{} {}_{}_{}\n{}\n".format(ppid, contig, start, end, seq[start:end]))
    out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Write a fasta file of prophage sequences from a set of genomes")
    parser.add_argument('-p', help='phispy directory', required=True)
    parser.add_argument('-g', help='genome directory (we look for the file "contigs")', required=True)
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    locs = parse_prophage_tbl(args.p)
    parse_contigs(locs, args.g, args.o)
