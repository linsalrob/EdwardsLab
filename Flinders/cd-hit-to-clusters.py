"""
Convert cdhit clusters + fasta file to a directory of fasta files with clusters
"""

import os
import sys
import argparse
import re
from roblib import read_fasta

__author__ = 'Rob Edwards'


def read_cd_hit_clusters(cdhitfile):
    """
    Read the clusters and return an array of arrays of clusters
    :param cdhitfile: the file of CD hit clusters to read
    :return: an array of arrays of clusers. Each element is the fasta id
    """

    results = []
    with open(cdhitfile, 'r') as f:
        cluster = []
        for l in f:
            if l.startswith('>'):
                if cluster:
                    results.append(cluster)
                cluster = []
            else:
                m = re.search('>(\d+)', l)
                if not m:
                    sys.stderr.write("No sequence id found in {}".format(l))
                    continue
                cluster.append(m.groups()[0])
                # sys.stderr.write(" ".join(cluster) + "\n")

    results.append(cluster)
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert cdhit + clusters to a directory of fasta files')
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-c', help='cd hit clusters file', required=True)
    parser.add_argument('-d', help='output directory', required=True)
    args = parser.parse_args()



    clusters = read_cd_hit_clusters(args.c)
    fa = read_fasta(args.f)
    if not os._exists(args.d):
        os.mkdir(args.d)

    for c in range(len(clusters)):
        with open(os.path.join(args.d, "Cluster{}.fasta".format(c)), 'w') as out:
            for seq in clusters[c]:
                out.write(">{}\n{}\n".format(seq, fa[seq]))

