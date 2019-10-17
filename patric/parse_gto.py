"""
Parse a GTO object
"""

import os
import sys
import argparse

from roblib import bcolors
import json


def list_keys(gto, verbose=False):
    """
    List the primary keys in the patric file
    :param gto: the json gto
    :param verbose: more output
    :return:
    """
    print("{}".format("\n".join(gto.keys())))

def dump_json(gto, k, verbose=False):
    """
    Print out the json representation of some data
    :param gto: the json gto
    :param k: the key to dump (none for everything)
    :param verbose: more output
    :return:
    """

    if k:
        if k in gto:
            print(json.dumps(gto, indent=4))
        else:
            sys.stderr.write(f"{bcolors.RED}ERROR: {k} not found.{bcolors.ENDC}\n")
    else:
        pp.pprint(f"{gto}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot a heatmap")
    parser.add_argument('-f', help='gto file', required=True)
    parser.add_argument('-l', help='list the primary keys and exit', action='store_true')
    parser.add_argument('-d', help='dump some part of the json object', action='store_true')
    parser.add_argument('-k', help='json primary key (e.g. for dumping, etc)')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    gto = json.load(open(args.f, 'r'))

    if args.l:
        list_keys(gto, args.v)
        sys.exit(0)

    if args.d:
        dump_json(gto, args.k, args.v)

