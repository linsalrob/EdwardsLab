"""
Parse a jplacer file and see what happens!
"""


import os
import sys
import argparse

import json
import re

from ete3 import Tree
from ete3.parser.newick import NewickError


def load_jplacer(jpf):
    """
    load the jplacer file and return the tree
    :param jpf: The jplacer file
    :return: the data structure of the tree
    """

    with open(jpf, 'r') as f:
        data = json.load(f)

    return data


def explore_jplacer(data):
    """
    Parse a jplacer data structure
    :param the data structure from the jplacer file:
    :return:
    """


    # print("{}\n".format(data.keys()))
    # sys.exit()
    print("{}".format(data['fields']))

    print("{}".format(data['placements']))
    sys.exit(0)

    """
    d1 = data['placements'][0]
    print("{}".format(d1))
    """

    for d in data['placements']:
        if 'nm' in d:
            print("{}".format(d))
            break
    """

    for d in data['placements']:
        print("\n".join(d.keys()))

    """

def get_placements(data, distmeasure):
    """
    Get the placements and return a dict with the keys being the edge numbers where to do
    the insertions and the values being a dict of nodes to insert at that point and their distances.

    For this purposes we multiple the distance measure by the modifier for each entry

    TODO: we have not implemented the approach for a single insertion as we don't have an example of that (yet!)

    :param data: the parsed jplacer tree
    :param distmeasure: the distance measure to use
    :return: a dict of placement edge_numbers and sets of ids to add
    """

    # first make sure the tree fields are in the correct order!
    posn = data['fields'].index('edge_num')
    if distmeasure not in data['fields']:
        sys.stderr.write("Crap. We do not have {} in our possible fields: {}\n".format(distmeasure, data['fields']))
        sys.exit(-1)
    distn = data['fields'].index(distmeasure)
    placements = {}

    for pl in data['placements']:
        multiplier = {}
        if 'n' in pl:
            sys.stderr.write("Crap, not sure what to do because I've never seen an example. You should be able to figure out from what I did with nm\n")
            sys.exit(-1)
        if 'nm' in pl:
            for i in pl['nm']:
                multiplier[i[0].replace(' ', '_')] = i[1]
        for p in pl['p']:
            edge_num = p[posn]
            distance = p[distn]
            if edge_num not in placements:
                placements[edge_num] = {}
            for thisid in multiplier:
                placements[edge_num][thisid] = multiplier[thisid] * distance * 1.0

    return placements


def parse_jplacer_tree(data):
    """
    Extract the tree from the jplacer data structure and make it into an ete3 object
    :param data: the jplacer data structure
    :return:
    """


    try:
        tree = Tree(data['tree'], quoted_node_names=True, format=1)
    except NewickError as n:
        tt = re.sub(r'(\:[\d\.]+){\d+}', r'\1', data['tree'])
        tt = re.sub(r'{\d+};$', ';', tt)
        tree = Tree(tt, quoted_node_names=True, format=1)

    return tree


def find_a_node(tree, nodeid):
    """
    Find a specific node in the tree
    :param tree: the tree to search
    :param nodeid: the node id to look for. This should be a string
    :type nodeid: str
    :return:
    """

    print("Traversing and looking for {{{}}}".format(nodeid))

    for t in tree.traverse("preorder"):
        if "{{{}}}".format(nodeid) in t.name:
            print("Found {}".format(t.name))
            t.add_child(name="NEW CHILDE")


    print(tree.write(format=1))

def insert_new_nodes(tree, placements, verbose=False):
    """
    Insert the new nodes in the tree at the correct place
    :param tree: The phylogenetic tree
    :param placements: the list of edges to place in the right place
    :return:
    """

    added = set()

    for t in tree.traverse("postorder"):
        m = re.search('{(\d+)}', t.name)
        if not m:
            continue
        thisid = int(m.groups()[0])

        if thisid in placements:
            added.add(thisid)
            for n in placements[thisid]:
                t.add_child(name=n, dist=placements[thisid][n])

    if verbose:
        for p in placements:
            if p not in added:
                sys.stderr.write("We did not find a node to add {}\n".format(p))

    return tree

def rename_leaves_taxids(tree):
    """
    Rename the leaf nodes with just the NCBI taxonomy ID if we have it
    :param tree: the tree to rename
    :return: the tree with renamed leaves
    """

    for n in tree.get_leaves():
        m = re.search(r'\[(\d+)\]', n.name)
        if m:
            n.name = m.groups()[0]
    return tree


def write_tree(tree, outputf):
    """
    Write the tree to a file.
    :param tree: The tree to write
    :param outputf: The output filename
    :return:
    """

    tree.write(outfile=outputf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse a jplacer file')
    parser.add_argument('-j', help='jplacer file', required=True)
    parser.add_argument('-o', help='output file to write the tree to', required=True)
    parser.add_argument('-d', help='Distance (can be either distal_length (default) or pendant_length', default='distal_length')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.d not in ['distal_length', 'pendant_length']:
        sys.stderr.write("sorry, at the moment the -d option must be either distal_length or pendant_length)\n")
        sys.exit()


    data = load_jplacer(args.j)

    tree = parse_jplacer_tree(data)

    # explore_jplacer(data)

    # find_a_node(tree, "2161")
    placements = get_placements(data, args.d)
    tree = insert_new_nodes(tree, placements, args.v)
    tree = rename_leaves_taxids(tree)

    write_tree(tree, args.o)
