"""
Just list clusters and number of members to enable sorting and summarizing!
"""

import os
import sys
import argparse
import re
__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def read_cd_hit_clusters(cdhitfile):
    """
    Read the clusters and return an array of arrays of clusters
    :param cdhitfile: the file of CD hit clusters to read
    :return: an array of arrays of clusers. Each element is the fasta id
    """

    results = {}
    lastid = None
    with open(cdhitfile, 'r') as f:
        for l in f:
            if l.startswith('>'):
                lastid = l.strip()
                results[lastid] = []
            else:
                m = re.search('>(\S+)', l)
                if not m:
                    sys.stderr.write("No sequence id found in {}".format(l))
                    continue
                r = m.groups()[0]
                while r.endswith('.'):
                    r = r[:-1]
                results[lastid].append(r)

    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-c', help='cd hit clusters file', required=True)
    args = parser.parse_args()

    c = read_cd_hit_clusters(args.c)
    for cl in c:
        print(f"{len(c[cl])}\t{cl}")