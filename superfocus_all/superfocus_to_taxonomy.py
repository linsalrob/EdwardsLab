"""
Parse the superfocus m8 diamon output files and sumamrize the taxonomies
"""

import os
import sys
import argparse
import re
from taxon import taxonomy_hierarchy_as_list, get_taxonomy_db
from roblib import colors

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

taxonomy = {}


def printmatches(seqid, matches):
    """
    Summarize and print he matches
    :param seqid: the sequence ID
    :param matches:  the matches array
    :return:
    """


    if len(matches) == 0:
        print(f"{seqid}\tUnknown")
        return

    ranks = ["Root", "", "", "", "", "", "", "", ""]
    for i in range(0, 8):
        t = set()
        for r in matches:
            t.add(r[i])
        if len(t) == 1:
            ranks[i + 1] = t.pop()
    s = "\t".join(ranks)
    print(f"{seqid}\t{s}")

def parse_m8(m8f, evalue, verbose=False):
    """
    Parse the m8 file and ... do something!
    :param m8f: the m8 output file from diamond
    :param evalue: the maximum evalue
    :param verbose: more output (maybe)
    :return:
    """

    fig = re.compile('fig\|(\d+)\.\d+')
    c = get_taxonomy_db()
    global taxonomy

    matches = []
    if verbose:
        sys.stderr.write(f"{colors.GREEN}Reading {m8f}{colors.ENDC}\n")
    lastid = None
    with open(m8f, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if float(p[10]) > evalue:
                continue
            if lastid and p[0] != lastid:
                printmatches(lastid, matches)
                matches = []
            lastid = p[0]
            m = fig.match(p[1])
            if m:
                tid = m.group(1)
                if tid not in taxonomy:
                    taxonomy[tid] = taxonomy_hierarchy_as_list(c, tid, True)
                matches.append(taxonomy[tid])


def parse_m8_tophit(m8f, evalue, verbose=False):
    """
    Parse the m8 file and ... do something!
    :param m8f: the m8 output file from diamond
    :param evalue: the maximum evalue
    :param verbose: more output (maybe)
    :return:
    """

    fig = re.compile('fig\|(\d+)\.\d+')
    c = get_taxonomy_db()
    global taxonomy

    matches = []
    if verbose:
        sys.stderr.write(f"{colors.GREEN}Reading {m8f}{colors.ENDC}\n")
    printed=set()
    with open(m8f, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if float(p[10]) > evalue:
                continue
            if p[0] in printed:
                continue
            m = fig.match(p[1])
            if m:
                tid = m.group(1)
                if tid not in taxonomy:
                    taxonomy[tid] = taxonomy_hierarchy_as_list(c, tid, True)
                r = "\t".join(taxonomy[tid])
                printed.add(p[0])
                print(f"{p[0]}\t{r}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-d', '--directory', help='Super focus output directory. Needs the m8 files', required=True)
    parser.add_argument('-e', '--evalue', help='max evalue. Default 1e-5', type=float, default=1e-5)
    parser.add_argument('-t', '--tophit', help='only use the tophit', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    for f in os.listdir(args.directory):
        if f.endswith('.m8'):
            if args.tophit:
                parse_m8_tophit(os.path.join(args.directory, f), args.evalue, args.v)
            else:
                parse_m8(os.path.join(args.directory, f), args.evalue, args.v)

