"""
Generate a color strip that includes Bacteria, Archaea, and Metagenomes
"""

import os
import sys
import argparse
import re
from roblib import stream_fastq
from taxon import get_taxonomy_db, get_taxonomy

c = get_taxonomy_db() 

def fq_ids(fnames, verbose=False):
    """
    Get a list of fastq ids for each of the files in fnames
    :param fnames: a list of files
    :return: a dict of ids
    """

    if verbose:
        sys.stderr.write("Reading fastq files\n")

    fqids = {}
    for f in fnames:
        for seqid, fullid, seq, qual in stream_fastq(f):
            # note we store several versions of the id as phylosift does some munging on them
            fqids[fullid] = f
            fullid = fullid.replace(' ', '_')
            fqids[fullid] = f

    return fqids


def determine_phylogeny(fn, fqids, verbose=False):
    """
    Determine if we know the phylogeny of this thing
    :param fn: the feature name
    :param fqids: the dict of ids->fastq files
    :return: the type (currently Bacteria, Archaea, Eukaryota, Metagenome, or unknown) and if a metagenome the type of metagenome
    """


    if fn in fqids:
        return "Metagenome", fqids[fn]
    m = re.sub('\.\d+\.\d+$', '', fn)
    if m in fqids:
        return "Metagenome", fqids[m]
    m = re.search('\[(\d+)\]', fn)
    if not m:
        if verbose:
            sys.stderr.write("There is no taxid in {} and it is not in the fastq file\n".format(fn))
        return "Unknown", None

    tid = m.groups()[0]
    t,n = get_taxonomy(tid, c)
    if not t:
        if verbose:
            sys.stderr.write("Can't find tax for {} in the db\n".format(tid))
        return "Unknown", None

    while t.parent > 1 and t.parent != 131567:
        # 131567 is cellular organisms
        t,n = get_taxonomy(t.parent, c)
    return n.scientific_name, None


def color_mg(leaff, fqfiles, inneroutputf, outeroutputf, readdeff, verbose=False):
    """
    Write the color file given the metagenome file and the dict of IDs
    :param leaff: the leaves file
    :param fqfiles: the list of fastq files
    :param readdeff: read definition file to write
    :param inneroutputf: the file to write the inner circle to
    :param outeroutputf: the file to write the outer circle to
    :return:
    """

    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']

    fqids = fq_ids(fqfiles, verbose)
    # get the list of everything
    inner = {}
    outer = {}
    domains = set()
    stypes = set()
    if verbose:
        sys.stderr.write("Determining phylogeny for all leaves\n")
    with open(leaff, 'r') as f:
        with open(readdeff, 'w') as readout:
            for l in f:
                l = l.strip()
                dom, nm = determine_phylogeny(l, fqids, verbose)
                if nm:
                    outer[l]=nm
                    stypes.add(nm)
                inner[l]=dom
                domains.add(dom)
                readout.write("{}\t{}\t{}\n".format(l, dom, nm))

    if len(domains) > len(colors):
        sys.stderr.write("Crap. We don't have enough colors. Inner has {}. Try http://colorbrewer2.org/ for more colors!\n".format(len(inner)))
        sys.exit(-1)

    if verbose:
        sys.stderr.write("We have {} domains\n".format(len(domains)))


    innername = list(domains)
    outername = list(stypes)

    with open(inneroutputf, 'w') as out:
        out.write("DATASET_COLORSTRIP\nSEPARATOR SPACE\n")
        out.write("DATASET_LABEL Kingdom\n")
        out.write("COLOR #ff0000\n")

        # write the legend
        out.write("LEGEND_TITLE Kingdom Legend\n")
        out.write("LEGEND_COLORS")
        for i, j in enumerate(innername):
            out.write(" {}".format(colors[i]))
        out.write("\n")
        out.write("LEGEND_SHAPES")
        for i, j in enumerate(innername):
            out.write(" 1")
        out.write("\n")
        out.write("LEGEND_LABELS")
        for i, j in enumerate(innername):
            out.write(" {}".format(j))
        out.write("\n")

        # place the strip
        out.write("STRIP_WIDTH 25\n")
        out.write("COLOR_BRANCHES 1\n")
        out.write("DATA\n")

        for l in inner:
            y = re.sub('[=\[\]:]', '_', l)
            c = colors[innername.index(inner[l])]
            out.write("{} {} {}\n".format(y, c, inner[l]))

    with open(outeroutputf, 'w') as out:
        out.write("DATASET_COLORSTRIP\nSEPARATOR SPACE\n")
        out.write("DATASET_LABEL Metagenome\n")
        out.write("COLOR #00ff00\n")

        # write the legend
        out.write("LEGEND_TITLE Metagenome Legend\n")
        out.write("LEGEND_COLORS")
        for i, j in enumerate(outername):
            out.write(" {}".format(colors[i]))
        out.write("\n")
        out.write("LEGEND_SHAPES")
        for i, j in enumerate(outername):
            out.write(" 2")
        out.write("\n")
        out.write("LEGEND_LABELS")
        for i, j in enumerate(outername):
            out.write(" {}".format(j))
        out.write("\n")

        # place the strip
        out.write("STRIP_WIDTH 25\n")
        out.write("COLOR_BRANCHES 1\n")
        out.write("DATA\n")

        for l in outer:
            y = re.sub('[=\[\]:]', '_', l)
            c = colors[outername.index(outer[l])]
            out.write("{} {} {}\n".format(y, c, outer[l]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Color a list of metagenomes')
    parser.add_argument('-l', help='leaves list file', required=True)
    parser.add_argument('-d', help='directory of fastq files')
    parser.add_argument('-f', help='fastq file(s) upon which to base the coloring', action='append')
    parser.add_argument('-i', help='inner circle output filename (the domains)', required=True)
    parser.add_argument('-o', help='outer circle output filename (the metagenomes)', required=True)
    parser.add_argument('-r', help='read domain output filename (to make downstream processing easier)', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    fqfiles = []
    if args.f:
        fqfiles = args.f
    if args.d:
        for q in os.listdir(args.d):
            if q.endswith('fastq'):
                fqfiles.append(os.path.join(args.d, q))
    if len(fqfiles) == 0:
        sys.stderr.write("You must supply some fastq files with either -d (directory) or -f (files)\n")
        sys.exit(-1)

    color_mg(args.l, fastqfiles, args.i, args.o, args.r, args.v)
