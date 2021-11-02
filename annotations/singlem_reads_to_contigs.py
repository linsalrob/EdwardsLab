"""
Merge the annotations from singlem to contigs.

Given:
    - singlem output
    - bam file with reads mapped to contigs

Return
    - a percentage annotation of each contig
    - perhaps also a definitive annotation
"""

import os
import sys
import argparse
import pysam

__author__ = 'Rob Edwards'


def read_singlem(singlemf, verbose=False):
    """
    Read the singlem file
    :param str singlemf: the singlem output file name
    :param bool verbose: more output
    :return dict: the mapping of reads -> taxonomy
    """
    if verbose:
        sys.stderr.write(f"Reading {singlemf}\n")

    singlem = {}
    with open(singlemf, 'r') as f:
        for l in f:
            if l.startswith('gene\t'):
                continue
            p = l.strip().split("\t")
            for rd in p[6].split():
                # strip off the read number because bowtie2 doesn't keep that
                if rd.endswith('/1') or rd.endswith('/2'):
                    rd = rd[:-2]
                if rd in singlem and singlem[rd] != p[5]:
                    sys.stderr.write(f"Read {rd} gets different taxonomy for /1 and /2\n")
                singlem[rd] = p[5]
    return singlem

def parse_bamfile(bamfile, mapping, verbose=False):
    """
    :param str bamfile: the name of the bamfile
    :param dict mapping: the mapping of reads to values (e.g. singlem, kraken, etc)
    :param bool verbose: more output

    """

    counts = {}
    sam = pysam.AlignmentFile(bamfile, "rb")
    for aln in sam.fetch():
        if not aln.is_secondary:
            read = aln.query_name
            contig = aln.reference_name
            mapq =  aln.mapping_quality # at the moment we don't use this but we could filter things
            if read not in mapping:
                # no point in generating an error because many reads will not be assigned
                continue
            if contig not in counts:
                counts[contig] = {}
            counts[contig][mapping[read]] = counts[contig].get(mapping[read], 0) + 1

    return counts

def resolve_counts(counts, outputfile, threshold=0.75, verbose=False):
    """
    Resolve and print the output

    :param counts: The dictionary of contigs and mapping and their counts
    :type counts: dict
    :param outputfile: the file base name to write
    :type outputfile: str
    :param float threshold: the cutoff for asserting this is the chosen object
    :param verbose: more output
    :type verbose: bool
    :return:
    :rtype:
    """

    if verbose:
        sys.stderr.write(f"Writing all matches to {outputfile}.allmatches.tsv  and best hits to {outputfile}.besthits.tsv")

    with open(f"{outputfile}.allmatches.tsv", 'w') as all, open(f"{outputfile}.besthits.tsv", 'w') as best:
        for contig in sorted(counts.keys()):
            total = sum(counts[contig].values())
            all.write(contig)
            taxas = []
            fracts = []
            for t in counts[contig]:
                all.write(f"\t{t}|{counts[contig][t]}")
                taxas.append(t)
                fracts.append(counts[contig][t] / total)
            all.write("\n")
            maxfr = max(fracts)
            if fracts.count(maxfr) > 1:
                if verbose:
                    sys.stderr.write(f"For {contig} two different data have the same maximum: {maxfr}\n")
                decision = 'underterminable'
            elif maxfr < threshold:
                if verbose:
                    sys.stderr.write(f"For {contig} the maximum fraction ({maxfr}) was below the threshold ({threshold})\n")
                decision = "Below threshold"
            else:
                decision = taxas[fracts.index(maxfr)]

            best.write(f"{contig}\t{decision}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract singlem annotations and apply them to contigs')
    parser.add_argument('-b', '--bam', required=True,
                        help='bam file with read names')
    parser.add_argument('-o', '--output', required=True,
                        help='basename for output. We add .allmatches.tsv and .besthits.tsv to this name')
    parser.add_argument('-s', '--singlem', required=True,
                        help='singlem output file that requires reads in column 5 (ie. you need the --output_extras flag')
    parser.add_argument('-t', '--threshold', default=0.75, type=float,
                        help='Threshold to decide this is the unique classification for the contig. Default [Default: %(default)d]')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    sm = read_singlem(args.singlem, args.verbose)
    counts = parse_bamfile(args.bam, sm, args.verbose)
    resolve_counts(counts, args.output, args.threshold, args.verbose)