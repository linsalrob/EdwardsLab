"""
Given a directory of fastq files, a classification file that describes what the fastq files are,
and a tree where the leaves have been placed into the tree (i.e. from parse_rename_write.py)
we make a file that has [id, domain, type, source]
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


def fq_classification(fqclass, verbose=False):
    """
    Read the fastq classification file
    :param fqclass: the classification file that has the file name and then arbitrary classifications separated by tabs
    :param verbose: more output
    :return: a dict of the classification. Guaranteed that all have the same number of elements.
    """

    classi = {}
    maxlen = 0
    with open(fqclass, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if len(p) > maxlen:
                maxlen = len(p)
            classi[p[0]] = p[1:]

    for i in classi:
        while len(classi[i]) < maxlen:
            classi[i].append("None")


    strclassi = {x:"\t".join(classi[x]) for x in classi}

    return strclassi


def determine_phylogeny(fn, fqids, verbose=False):
    """
    Determine if we know the phylogeny of this thing
    :param fn: the feature name
    :param fqids: the dict of ids->fastq files
    :return: the type (currently Bacteria, Archaea, Eukaryota, Metagenome, or unknown), the id in the fastq file, and if a metagenome the type of metagenome
    """


    if fn in fqids:
        return "Metagenome", fn, fqids[fn]
    m = re.sub('\.\d+\.\d+$', '', fn)
    if m in fqids:
        return "Metagenome", m, fqids[m]
    m = re.sub('\.[\d\.]+$', '', fn)
    if m in fqids:
        return "Metagenome", m, fqids[m]
    m = re.search('\[(\d+)\]', fn)
    if not m:
        if verbose:
            sys.stderr.write("There is no taxid in {} and it is not in the fastq file\n".format(fn))
        return "Unknown", fn, None

    tid = m.groups()[0]
    t,n = get_taxonomy(tid, c)
    if not t:
        if verbose:
            sys.stderr.write("Can't find tax for {} in the db\n".format(tid))
        return "Unknown", fn, None

    while t.parent > 1 and t.parent != 131567:
        # 131567 is cellular organisms
        t,n = get_taxonomy(t.parent, c)
    return n.scientific_name, fn, None


def write_output(leaff, fqfiles, classifile, readdeff, verbose=False):
    """
    Write an output file that categorizes each leaf
    :param leaff: the leaves file
    :param fqfiles: the list of fastq files
    :param classifile: the classification file
    :param readdeff: read definition file to write
    :return:
    """

    cl = fq_classification(classifile, verbose)
    fqids = fq_ids(fqfiles, verbose)
    # get the list of everything
    domains = set()
    stypes = set()
    if verbose:
        sys.stderr.write("Determining phylogeny for all leaves\n")
    with open(leaff, 'r') as f:
        with open(readdeff, 'w') as readout:
            for l in f:
                l = l.strip()
                dom, id_in_fq, nm = determine_phylogeny(l, fqids, verbose)
                if nm:
                    # this also means that l is in fqids, so we can get the classification
                    thisfq = fqids[id_in_fq].split(os.path.sep)[-1]
                    print("Looking for thisfq: {}".format(thisfq))
                    clstr=""
                    if thisfq not in cl:
                        sys.stderr.write(f"ERROR: {thisfq} not found in the fastq classification file\n")
                    else:
                        clstr = cl[thisfq]
                    readout.write("{}\t{}\t{}\t{}\t{}\n".format(l, id_in_fq, dom, nm, clstr))
                else:
                    readout.write("{}\t{}\t{}\n".format(l, l, dom))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Color a list of metagenomes')
    parser.add_argument('-l', help='leaves list file', required=True)
    parser.add_argument('-c', help='fastq classification file', required=True)
    parser.add_argument('-d', help='directory of fastq files')
    parser.add_argument('-f', help='fastq file(s) [one or more can be specified]', action='append')
    parser.add_argument('-o', help='output file to write to', required=True)
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

    write_output(args.l, fqfiles, args.c, args.o, args.v)
