"""
Just list the products of all the features in a genbank file
"""

import os
import sys
import argparse

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

    record = SeqioFilter(SeqIO.parse(gbf, "genbank"))
    for entry in record:
        for feature in entry.get_features('CDS'):
            print(f"{entry.id}\t{feature.id}\t{feature.function}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='genbank file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    list_features(args.f, args.v)
