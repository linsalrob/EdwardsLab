"""
Figure out if proteins are phage, hypothetical, etc
"""

import os
import sys
import argparse
from PhiSpyModules import is_phage_func, is_unknown_func, is_gzip_file
from roblib import is_hypothetical, seqio_filter
import gzip
from Bio import Seq, SeqIO

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

def list_features(gbf, verbose=False):
    """
    List all the features in gbf
    :param gbf: genbank file
    :param verbose: more output
    :return:
    """

    if verbose:
        print(f"Parsing {gbf}", file=sys.stderr)

    try:
        if is_gzip_file(gbf):
            handle = gzip.open(gbf, 'rt')
        else:
            handle = open(gbf, 'r')
    except IOError as e:
        print(f"There was an error reading {args_parser.infile}: {e}", file=sys.stderr)
        sys.exit(20)

    record = seqio_filter.SeqioFilter(SeqIO.parse(handle, "genbank"))
    for entry in record:
        hypo = 0
        phage = 0
        notp  = 0
        for feature in entry.get_features('CDS'):
            fn = feature.function
            if is_unknown_func(fn) or is_hypothetical(fn):
                hypo += 1
            elif is_phage_func(fn):
                phage += 1
            else:
                notp += 1

        total = hypo + phage + notp
        print(f"{entry.id}\t{phage:,}\t{hypo:,}\t{notp:,}\t{total:,}")

        handle.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Categorize phage genes')
    parser.add_argument('-g', help='genbank input file', required=True, action='append')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    print("File\tPhage\tHypothetical\tNot Phage\tTotal")
    for g in args.g:
        list_features(g, args.v)

