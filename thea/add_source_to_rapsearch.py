"""

"""

import os
import sys
import argparse
import gzip

__author__ = 'Rob Edwards'

def read_sources(sf, verbose=False):
    """
    Read the sources.txt file
    :param sf: read the sources
    :param verbose:
    :return:
    """

    source={}
    with open(sf, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            source[p[0]]=p[1]
    return source

def read_rapsearch_out(source, rpfile, ofile, verbose=False):
    """
    Read the rapsearch output file
    :param source:
    :param rpfile:
    :param verbose:
    :return:
    """

    if rpfile.endswith('.gz'):
        qin = gzip.open(rpfile, 'rt')
    else:
        qin = open(rpfile, 'r')

    with open(ofile, 'w') as out:
        for l in qin:
            if l.startswith("#"):
                out.write(l)
                continue
            p=l.strip().split("\t")
            if p[1] not in source:
                sys.stderr.write("No source for {}\n".format(p[1]))
                out.write("\t{}".format(l))
                continue
            out.write("{}\t{}".format(source[p[1]], l))
    qin.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add the source of the file to rapsearch')
    parser.add_argument('-f', help='input rapsearch m8 file', required=True)
    parser.add_argument('-o', help='output file of rapsearch m8 with additional column', required=True)
    parser.add_argument('-s', help='sources list. Default=sources.txt', default="sources.txt")
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    s = read_sources(args.s, args.v)
    read_rapsearch_out(s, args.f, args.o, args.v)