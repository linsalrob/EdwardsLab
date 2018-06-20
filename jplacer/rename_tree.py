"""

Rename the internal nodes in the pplacer tree using my awesome NCBI taxonomy
"""

import os
import sys
import argparse
import json
import re

from ete3 import Tree
from ete3.parser.newick import NewickError

from taxon import get_taxonomy_db, get_taxonomy

def load_jplacer(jpf, verbose=False):
    """
    load the jplacer file and return the tree
    :param jpf: The jplacer file
    :return: the data structure of the tree
    """

    with open(jpf, 'r') as f:
        data = json.load(f)

    return data


def parse_jplacer_tree(data, verbose=False):
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

def rename_nodes(tree, verbose=False):
    """
    Rename the nodes based on everything below me
    """

    # connect to the SQL dataabase
    c = get_taxonomy_db()

    wanted_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
    wanted_levels.reverse() # too lazy to write in reverse :)
    taxonomy = {}
    # first get all the leaves and their parents. This is just to speed things up ... maybe
    for l in tree.get_leaves():
        m = re.search('\[(\d+)\]', l.name)
        if not m:
            if verbose:
                sys.stderr.write("No taxid in {}\n".format(l.name))
            continue
        tid = m.groups()[0]
        taxonomy[l.name] = {}
        t,n = get_taxonomy(tid, c)
        if not t:
            continue
        while t.parent != 1 and t.taxid != 1:
            if t.rank in wanted_levels:
                taxonomy[l.name][t.rank] = n.scientific_name
            t,n = get_taxonomy(t.parent, c)
        
    # now traverse every node that is not a leaf and see if we can some up with a 
    # unique name for the node!
    sys.stderr.write("Traversing the tree\n")
    for n in tree.traverse("preorder"):
        if n.is_leaf():
            continue
        sys.stderr.write("Checking {}\n".format(n.name))
        taxs = {w:set() for w in wanted_levels}
        for l in n.get_leaves():
            if l.name not in taxonomy:
                continue
            for w in wanted_levels:
                if w in taxonomy[l.name]:
                    taxs[w].add(taxonomy[l.name][w])
        # which is the LOWEST level with a single taxonomy
        for w in wanted_levels:
            if len(taxs[w]) == 1:
                newname = "{} (rank: {})".format(taxs[w].pop(), w)
                if verbose:
                    True
                sys.stderr.write("Changing name from: {} to {}\n".format(n.name, newname))
                n.name = newname
                break
    return tree

        
    






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse a jplacer file')
    parser.add_argument('-j', help='jplacer file', required=True)
    parser.add_argument('-o', help='output file to write the tree to', required=True)
    parser.add_argument('-v', help='verbose', action='store_true')
    args = parser.parse_args()

    data = load_jplacer(args.j, args.v)

    tree = parse_jplacer_tree(data, args.v)

    tree = rename_nodes(tree, args.v)



