"""
Read the cophenetic matrix and add some additional metadata to it.
"""

import os
import sys
import argparse
import re

def read_domains(df):
    """
    Read and return the domains
    :param df: the domains file to read
    :return: the domain for each of the species

    """

    data = {}
    with open(df, 'r') as f:
        for l in f:
            l = re.sub('sharks_cat.fastq', 'Sharks', l)
            l = re.sub('fish_cat.fastq', 'Fish', l)
            l = re.sub('ts_cat', 'Thresher_shark', l)
            l = re.sub('ws_cat', 'Whale_shark', l)
            l = re.sub('ls_cat', 'Leopard_shark', l)
            l = re.sub('_cat', '', l)
            p = l.strip().split("\t")
            p[0] = re.sub('[=\[\]:]', '_', p[0])
            data[p[0]] = [p[1], p[2], None]

    return data

def read_ann(af, data):
    """
    Read the annotationsand return that as a hash
    :param af: anntotations file
    :param data: parsed data from domains files
    :return:
    """

    with open(af, 'r') as f:
        for l in f:
            l = re.sub('sharks_cat.fastq', 'Sharks', l)
            l = re.sub('fish_cat.fastq', 'Fish', l)
            l = re.sub('ts_cat', 'Thresher_shark', l)
            l = re.sub('ws_cat', 'Whale_shark', l)
            l = re.sub('ls_cat', 'Leopard_shark', l)
            l = re.sub('_cat', '', l)
            m = re.search('^(\S)+\s.*\s(\S+)$', l)
            if m:
                if m.groups()[0] in data:
                    data[groups()[0]][2] = m.groups()[1]
                else:
                    sys.stderr.write("{} not found in data\n".format(m.groups()[0]))
            else:
                sys.stderr.write("Can't parse {} from {}\n".format(l, af))
    return data

def add_to_cp_matrix(cpm, data, outputf):
    """
    write the data to the cophenetic matrix
    :param cpm: the cp matrix file
    :param data: the data to add
    :param outputf: the output file to write
    :return:
    """

    with open(cpm, 'r') as f:
        headerl = f.readline()
        with open(outputf, 'w') as out:
            out.write("Domain\tType\tSpecies\tSeq ID\t{}".format(headerl))
            for l in f:
                p = l.strip().split("\t")
                if p[0] in data:
                    out.write("{}\t{}".format("\t".join(data[p[0]]), l))
                else:
                    sys.stderr.write("No data for {}\n".format(p[0]))
                    out.write("\t\t\t{}".format(l))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Add data to the cophenetic matrix")
    parser.add_argument('-m', help='cophenetic matrix', required=True)
    parser.add_argument('-r', help='reads domains output (the -r output from generate_color_strip.py)', required=True)
    parser.add_argument('-a', help='annotations file (the -o file from fastq2color_strip.py)', required=True)
    parser.add_argument('-o', help='output filename', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    data = read_domains(args.r)
    data = read_ann(args.a, data)
    add_to_cp_matrix(args.m, data, args.o)