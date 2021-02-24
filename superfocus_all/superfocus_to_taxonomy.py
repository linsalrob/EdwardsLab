"""
Parse the superfocus m8 diamond output files and sumamrize the taxonomies
"""

import os
import sys
import argparse
import re
from taxon import taxonomy_hierarchy_as_list, get_taxonomy_db
from taxon.Error import EntryNotInDatabaseError
from roblib import colors

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

taxonomy = {}

ignore_tids = {'6666666', '88888881'}


def printmatches(seqid, matches):
    """
    Summarize and print the matches
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

def parse_m8(m8f, evalue, database=None, verbose=False):
    """
    Parse the m8 file and ... do something!
    :param m8f: the m8 output file from diamond
    :param evalue: the maximum evalue
    :param verbose: more output (maybe)
    :return:
    """

    fig = re.compile('fig\|(\d+)\.\d+')
    if database:
        c = get_taxonomy_db(database)
    else:
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
                if tid in ignore_tids:
                    continue
                if tid not in taxonomy:
                    try:
                        taxonomy[tid] = taxonomy_hierarchy_as_list(c, tid, True)
                    except EntryNotInDatabaseError as e:
                        ignore_tids.add(tid)
                        continue
                if taxonomy[tid]:
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
                if tid in ignore_tids:
                    continue
                if tid not in taxonomy:
                    try:
                        taxonomy[tid] = taxonomy_hierarchy_as_list(c, tid, True)
                    except EntryNotInDatabaseError as e:
                        ignore_tids.add(tid)
                        continue
                if taxonomy[tid]:
                    if 'metagenome' in taxonomy[tid][6].lower():
                        taxonomy.pop(tid)
                        ignore_tids.add(tid)
                        continue
                    r = "\t".join(taxonomy[tid])
                    printed.add(p[0])
                    print(f"{p[0]}\t{r}\ttaxid: {tid}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', '--file', help='Super focus m8 output file', required=True)
    parser.add_argument('-e', '--evalue', help='max evalue. Default 1e-5', type=float, default=1e-5)
    parser.add_argument('-t', '--tophit', help='only use the tophit', action='store_true')
    parser.add_arguemtn('-d', '--db', help='SQLite3 taxonomy database location')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    if args.tophit:
        parse_m8_tophit(args.file, args.evalue, args.v)
    else:
        parse_m8(args.file, args.evalue, args.db, args.v)

