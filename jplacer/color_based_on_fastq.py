"""
Color a list of metagenomes based on fastq files
"""

import os
import sys
import argparse
import re
from roblib import stream_fastq

def fq_ids(fnames):
    """
    Get a list of fastq ids for each of the files in fnames
    :param fnames: a list of files
    :return: a dict of ids
    """

    res = {}
    for f in fnames:
        for seqid, fullid, seq, qual in stream_fastq(f):
            # note we store several versions of the id as phylosift does some munging on them
            res[fullid] = f
            fullid = fullid.replace(' ', '_')
            res[fullid] = f

    return res


def color_mg(mgf, fnames, outputf):
    """
    Write the color file given the metagenome file and the dict of IDs
    :param mgf: the metagenome file
    :param fnames: the list of fastq files
    :param outputf: the file to write to
    :return:
    """

    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']
    if len(fnames) > len(colors):
        sys.stderr.write("Crap. We don't have enough colors. Try http://colorbrewer2.org/ for more colors!\n")
        sys.exit(-1)

    fids = fq_ids(fnames)


    with open(mgf, 'r') as f:
        with open(outputf, 'w') as out:
            out.write("TREE_COLORS\nSEPARATOR SPACE\nDATA\n")
            for l in f:
                l = l.strip()
                if l in fids:
                    y = l.replace(':', '_')
                    c = colors[fnames.index(fids[l])]
                    out.write("{} clade {} normal 10\n".format(y, c))
                else:
                    m = re.sub('\.\d+\.\d+$', '', l)
                    if m in fids:
                        y = l.replace(':', '_')
                        c = colors[fnames.index(fids[m])]
                        out.write("{} clade {} normal 10\n".format(y, c))
                    else:
                        sys.stderr.write("Error |{}| or |{}| not found in fastq ids\n".format(l, m))
                        continue






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Color a list of metagenomes')
    parser.add_argument('-m', help='metagenome list file', required=True)
    parser.add_argument('-f', help='fastq file(s) upon which to base the coloring', action='append', required=True)
    parser.add_argument('-o', help='output filename', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    color_mg(args.m, args.f, args.o)

