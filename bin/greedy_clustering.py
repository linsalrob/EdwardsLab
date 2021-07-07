"""
Given the output of correlations.py that calculates pairwise pearson correlations, generate all
possible clusters in a leader clustering approach
"""

import os
import sys
import argparse
from typing import List, Any, Dict, Set


def merge_clusters(clx, mmx, x, y, verbose=False):
    """
    Merge the clusters
    :param verbose: more output
    :type verbose: bool
    :param y: second member
    :type y: str
    :param x: first member
    :type x: str
    :param mmx: For a cluster ID, a list of the members
    :type mmx: Dict[int, Set[str]]
    :param clx: For the member of a cluster, the cluster number
    :type clx: dict[str, int]
    :rtype: dict[str, int], Dict[int, Set[str]]

    """
    if verbose:
        sys.stderr.write(f"Merging clusters with {x} and {y}\n")

    xc = clx[x]
    yc = clx[y]

    if xc == yc:
        # nothing to do!
        return clx, mmx

    for m in mmx[yc]:
        clx[m] = xc
    del mmx[yc]

    return clx, mmx


def cluster(inputfile, threshold, verbose=False):
    """
    Calculate the clusters
    :param str inputfile: the input file
    :param float threshold: the threshold for being included in a cluster
    :param verbose:  more output
    :return: the clusters and their members
    :rtype: dict[str, int], Dict[int, Set[str]]
    """

    clusters: Dict[str, int] = {}
    members: Dict[int, Set[str]] = {}
    cluster_number: int = 0

    with open(inputfile, 'r') as f:
        for li in f:
            p = li.strip().split("\t")
            if float(p[2]) < threshold:
                continue
            # is either p[0] or p[1] in a cluster?
            if p[0] in clusters and p[1] in clusters:
                if clusters[p[0]] != clusters[p[1]]:
                    clusters, members = merge_clusters(clusters, members, p[0], p[1], verbose)
            elif p[0] in clusters:
                xc = clusters[p[0]]
                clusters[p[1]] = xc
                members[xc].add(p[1])
            elif p[1] in clusters:
                xc = clusters[p[1]]
                clusters[p[0]] = xc
                members[xc].add(p[0])
            else:
                cluster_number += 1
                clusters[p[0]] = cluster_number
                clusters[p[1]] = cluster_number
                members[cluster_number] = {p[0], p[1]}

    return clusters, members


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', '--file', help='input file', required=True)
    parser.add_argument('-t', '--threshold', help='threshold', type=float, required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    clusters, members = cluster(args.file, args.threshold, args.verbose)
    with open(args.output, 'w') as out:
        for m in members:
            mems = "\t".join(members[m])
            out.write(f"{m}\t{mems}\n")
