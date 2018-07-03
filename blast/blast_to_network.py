"""
Create a networkx network from tab separated blast output.

Currently the default is for the edge weight to be the percent id * alignment length
"""

import os
import sys
import argparse
import networkx as nx

def load_blast(blastf, maxeval=None, minper=None, minlen=None):
    """
    Load a blast file into the network
    :param blastf: blast tab separated file to read
    :param maxeval: maximum E value to include in the network
    :param minper: minimum percent to include in the network
    :param minlen: minimum alignment length to include in the network
    :return: a networkx objext
    :return type: networkx
    """

    G = nx.Graph()
    with open(blastf, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            if maxeval and float(p[10]) > maxeval:
                continue
            if minper and float(p[2]) < minper:
                continue
            if minlen and float(p[3]) < minlen:
                continue

            weight = (float(p[2])/100) * float(p[3])
            G.add_edge(p[0], p[1], {'weight' : weight})

    return G



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create a network from blast tab separated output")
    parser.add_argument('-b', help='blast output file (required)', required=True)
    parser.add_argument('-e', help='maximum e-value', type=float)
    parser.add_argument('-p', help='minimum percent identity', type=float)
    parser.add_argument('-l', help='minimum alignment length', type=int)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()


    net = load_blast(args.b, maxeval=args.e, minper=args.p, minlen=args.l)

    print("Network:\n\tEdges: {}\n\tNodes: {}\n")