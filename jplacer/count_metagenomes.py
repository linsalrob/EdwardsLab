"""
Count metagenomes under phylum level
"""

import os
import sys
import argparse
from ete3 import Tree


def read_mg_types(mgfile):
    """
    Read the metagenome types from a file where the first four columns
    are [domain, type, species, seqid]. For example, sharks_fish.distance.labelled.tsv
    which is the output from join.pl
    :param mgfile: the sharks_fish.distance.labelled.tsv file
    :return: a dict of all the metagenomes
    """
    
    val = {}
    with open(mgfile, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            if "Metagenome" in p[0]:
                val[p[3]] = [p[1], p[2]]
    return val


def count_mg(treefile, mgdesc):
    """
    Read the tree and count below phylum
    :param treefile: the file with the tree
    :param mgdesc: the description of the metagenomes from read_mg_types
    """

    tree = Tree(treefile, quoted_node_names=True, format=1)

    ### NOTE TO SELF:
    """
    Probably want to make a multibar chart here, of the counts of fish/shark or 
    different types of fish/shark in our datasets.

    This part will so far get the counts. Just need to know how to print them
    """


    for n in tree.traverse("preorder"):
        if "r_phylum" in n.name:
            fishcount = {}
            typecount = {}
            for l in n.get_leaves():
                if l.name in mgdesc:
                    vals = mgdesc[l.name]
                    fishcount[vals[0]] = fishcount.get(vals[0], 0) + 1
                    typecount[vals[1]] = typecount.get(vals[1], 0) + 1
            for f in fishcount:
                print("{}\t{}\t{}".format(n.name, f, fishcount[f]))
            for t in typecount:
                print("{}\t{}\t{}".format(n.name, t, typecount[t]))




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count metagenomes on the tree')
    parser.add_argument('-t', help='tree file', required=True)
    parser.add_argument('-d', help='metagenome description file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    mt = read_mg_types(args.d)
    count_mg(args.t, mt)
