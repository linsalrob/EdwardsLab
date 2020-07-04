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

def parse_m8(m8f, eval, verbose=False):
    """
    Parse the m8 file and ... do something!
    :param m8f: the m8 output file from diamond
    :param eval: the maximum evalue
    :param verbose: more output (maybe)
    :return:
    """

    fig = re.compile('fig\|(\d+)\.\d+')
    c = get_taxonomy_db()
    global taxonomy

    matches = {}
    if verbose:
        sys.stderr.write(f"{colors.GREEN}Reading {m8f}{colors.ENDC}\n")
    with open(m8f, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if float(p[10] > eval):
                continue
            m = fig.match(p[1])
            matches[p[0]] = []
            if m:
                tid = m.group(1)
                if tid not in taxonomy:
                    taxonomy[tid] = taxonomy_hierarchy_as_list(c, tid, True)
                matches[p[0]].append(taxonomy[tid])

    for m in matches:
        if len(matches[m]) == 0:
            print(f"{m}\tUnknown")
            continue

        ranks = ["Root", "", "", "", "", "", "", "", ""]
        for i in range(0,8):
            t = set()
            for r in matches[m]:
                t.add(r[i])
            if len(t) == 1:
                ranks[i+1] = t.pop()
        s = "\t".join(ranks)
        print(f"{m}\t{s}")








if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-d', '--database', help='Super focus output directory. Needs the m8 files', required=True)
    parser.add_argument('-e', '--evalue', help='max evalue. Default 1e-5', type=float, default=1e-5)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    for f in os.listdir(args.d):
        if f.endswith('.m8'):
            parse_m8(os.path.join(args.d, f), args.e, args.v)

