"""
For an m8 blast file, read the second column as the genome id from PATRIC (fig), convert that
to a taxonomy, and then get the taxonomic profile of that sequence.

For all matches to a sequence that are better than our e-val, figure out the taxonomy that
covers all the matches.
"""

import os
import sys
import argparse
import re
from taxon import get_taxonomy_db, get_taxonomy


def id_from_blastfile(blastfile, evalue, verbose=False):
    """
    Read the blast file and yield a single sequence as a query and the match to its taxonomy id
    :param blastfile: the blast file to read
    :param evalue: the minimum evalue (col 10)
    :param verbose: more output
    :return: yields a query ID and a taxonomy ID
    """

    f = open(blastfile, 'r')

    while True:
        l = f.readline()
        if not l:
            break
        p = l.strip().split("\t")
        if float(p[10]) > evalue:
            continue
        m = re.search('fig\|(\d+)\.\d+', p[1])
        if not m:
            sys.stderr.write("Can't parse a taxonomy from {}\n".format(p[1]))
            continue
        tid = m.groups()[0]
        if tid.startswith('666.'):
            if verbose:
                sys.stderr.write("Skipped {}\n".format(tid))
            continue
        yield p[0], tid



def tid_to_tax_set(tid, verbose=False):
    """
    Convert a taxonomy ID to a set that contains the hierarchy.
    :param tid: the taxonomy ID
    :param verbose: print more stuff
    :return: the set of taxonomies
    """

    # connect to the SQL dataabase
    c = get_taxonomy_db()

    wanted_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
    taxonomy = set()
    t, n = get_taxonomy(tid, c)
    if not t:
        if verbose:
            sys.stderr.write("No taxonomy for {}\n".format(tid))
        return None

    while t.parent != 1 and t.taxid != 1:
        if t.rank in wanted_levels:
            taxonomy.add(n.scientific_name + ":" + t.rank)
        t, n = get_taxonomy(t.parent, c)

    return taxonomy


def taxa_sets(blastf, eval, verbose):
    """
    Define the taxonomy sets for a blast file
    :param blastf: the blast file to parse
    :param eval: the maximum evalue
    :param verbose: print more stuff
    :return: prints the id and set each time a new id is found
    """

    taxset = set()
    lastquery = None
    for query, tid in id_from_blastfile(blastfile=blastf, evalue=eval, verbose=verbose):
        ts = tid_to_tax_set(tid=tid, verbose=verbose)
        if None == ts:
            if verbose:
                sys.stderr.write("No known tax id for {} from {}. Skipped\n".format(query, tid))
            continue

        if query != lastquery:
            if lastquery and taxset:
                wanted_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
                taxalist = wanted_levels
                for i, w in enumerate(wanted_levels):
                    for t in taxset:
                        if f":{w}" in t:
                            taxalist[i] = t.replace(f":{w}", "")

                print("\t".join([lastquery, "; ".join(taxalist)]))
            taxset = ts
            lastquery = query
            continue

        taxset = taxset & ts



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Do some taxonomy-fu on blast from patric")
    parser.add_argument('-b', help='blast output file', required=True)
    parser.add_argument('-e', help='E value threshold default 1e-5', type=float, default=1e-5)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()


    taxa_sets(args.b, args.e, args.v)