"""
Explore the jplacer tree
"""

import os
import sys
import argparse
from ete3 import Tree
import re

def load_tree(treefile):
    """
    Just read the tree file
    :param treefile:
    :return:
    """

    return Tree(treefile, quoted_node_names=True, format=1)


def print_taxids(tree):
    """
    Print all internal and external nodes
    :param tree:
    :return:
    """

    for n in tree.traverse("postorder"):
        m = re.search('\[(\d+)\]', n.name)
        if not m:
            sys.stderr.write("No taxid in {}\n".format(n.name))
        else:
            print("{}\t{}".format(m.groups()[0], n.name))


def explore_tree(tree):
    """
    Explore a tree file in some way
    :param treefile: The tree to explore
    :return:
    """

    salmonellaleaf = None
    print("LEAVES:\n")
    for l in tree.get_leaves():
        print(l.name)
        if "497977" in l.name:
            salmonellaleaf = l

    print('\n\n')

    # start with one leave and explore the tree
    print(salmonellaleaf.name)
    for a in salmonellaleaf.get_ancestors():
        print(a.name)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a tree into a distance matrix')
    parser.add_argument('-t', help='Tree file', required=True)
    args = parser.parse_args()

    tree = load_tree(args.t)
    print_taxids(tree)