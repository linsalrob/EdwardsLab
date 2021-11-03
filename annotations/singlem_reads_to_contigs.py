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


def read_kraken(krakenf, verbose=False):
    """
    Parse the kraken output file.
    :param str krakenf: the kraken output file
    :param bool verbose:
    :return dict: the kraken reads and counts
    """

    if args.verbose:
        sys.stderr.write(f"Parsing Kraken file {krakenf}\n")

    kraken = {}
    with open(krakenf, 'r') as f:
        for l in f:
            if l.startswith('U'):
                # unclassified read
                continue
            p = l.strip().split("\t")
            if p[0].endswith('/1') or p[0].endswith('/2'):
                p[0] = p[0][:-2]
            if p[0] in kraken:
                continue
            kraken[p[0]] = p[2]
    return kraken

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
                    # sys.stderr.write(f"Read {rd} gets different taxonomy for /1 and /2\n")
                    continue
                singlem[rd] = p[5]
    return singlem

def parse_bamfile(bamfile, mapping, verbose=False):
    """
    :param str bamfile: the name of the bamfile
    :param dict mapping: a dict of methods that are dicts the mapping of reads to values (e.g. singlem, kraken, etc)
    :param bool verbose: more output

    """

    if verbose:
        sys.stderr.write(f"Parsing bam file {bamfile}\n")

    counts = {}
    for mthd in mapping:
        counts[mthd] = {}

    sam = pysam.AlignmentFile(bamfile, "rb")
    for aln in sam.fetch():
        if not aln.is_secondary:
            read = aln.query_name
            contig = aln.reference_name
            mapq =  aln.mapping_quality # at the moment we don't use this but we could filter things

            for mthd in mapping:
                if read in mapping[mthd]:
                    if contig not in counts[mthd]:
                        counts[mthd][contig] = {}
                    counts[mthd][contig][mapping[mthd][read]] = counts[mthd][contig].get(mapping[mthd][read], 0) + 1

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

    all_contigs = {}
    for mthd in counts:
        with open(f"{outputfile}.{mthd}.allmatches.tsv", 'w') as all, open(f"{outputfile}.{mthd}.besthits.tsv", 'w') as best:
            for contig in sorted(counts[mthd].keys()):
                total = sum(counts[mthd][contig].values())
                all.write(contig)
                taxas = []
                fracts = []
                for t in counts[mthd][contig]:
                    all.write(f"\t{t}|{counts[mthd][contig][t]}")
                    taxas.append(t)
                    fracts.append(counts[mthd][contig][t] / total)
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

                if contig not in all_contigs:
                    all_contigs[contig] = {}
                all_contigs[contig][mthd] = decision

                best.write(f"{contig}\t{decision}\n")

    with open(f"{outputfile}.comparison.tsv", 'w') as all:
        methods = sorted(counts.keys())
        all.write("\t".join(["contig"] + methods))
        all.write("\n")
        for c in all_contigs:
            all.write(c)
            for mthd in methods:
                all.write(f"\t{all_contigs[c].get(mthd, '')}")
            all.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract singlem annotations and apply them to contigs')
    parser.add_argument('-b', '--bam', required=True,
                        help='bam file with read names')
    parser.add_argument('-o', '--output', required=True,
                        help='basename for output. We add .allmatches.tsv and .besthits.tsv to this name')
    parser.add_argument('-s', '--singlem', required=True,
                        help='singlem output file that requires reads in column 5 (ie. you need the --output_extras flag')
    parser.add_argument('-k', '--kraken', help='kraken output file. The output from --output flag of kraken')
    parser.add_argument('-t', '--threshold', default=0.75, type=float,
                        help='Threshold to decide this is the unique classification for the contig. Default [Default: %(default)d]')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    data = {}
    if args.singlem:
        data['singlem'] = read_singlem(args.singlem, args.verbose)
    if args.kraken:
        data['kraken'] = read_kraken(args.kraken, args.verbose)
    if not data:
        sys.stderr.write("Please provide at least one or more classification files (kraken, singlem, etc)\n")
        sys.exit(0)
    counts = parse_bamfile(args.bam, data, args.verbose)
    resolve_counts(counts, args.output, args.threshold, args.verbose)