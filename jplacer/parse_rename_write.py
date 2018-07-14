"""
Parse a jplacer file, move our nodes into the tree, rename the branches, and then write a tree
"""


import os
import sys
import argparse

import json
import re

from ete3 import Tree
from ete3.parser.newick import NewickError

from taxon import get_taxonomy_db, get_taxonomy

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
def clean_newick_id(name):
    """
    Return a version of name suitable for placement in a newick file
    :param name: The name to clean up
    :return: a name with no colons, spaces, etc
    """
    name = name.replace(' ', '_')
    name = name.replace(':', '_')

    return name

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
                newid = clean_newick_id(thisid)
                placements[edge_num][newid] = multiplier[thisid] * distance * 1.0

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
            t.add_child(name="NEW CHILD")


    print(tree.write(format=1))

def insert_new_nodes(tree, placements, namesf=None, verbose=False):
    """
    Insert the new nodes in the tree at the correct place
    :param tree: The phylogenetic tree
    :param placements: the list of edges to place in the right place
    :return:
    """

    added = set()
    addednames = set()


    for t in tree.traverse("postorder"):
        m = re.search('{(\d+)}', t.name)
        if not m:
            continue
        thisid = int(m.groups()[0])

        if thisid in placements:
            added.add(thisid)
            for n in placements[thisid]:
                addednames.add(n)
                t.add_child(name=n, dist=placements[thisid][n])

    if verbose:
        for p in placements:
            if p not in added:
                sys.stderr.write("We did not find a node to add {}\n".format(p))
    if namesf:
        with open(namesf, 'w') as namesout:
            for p in addednames:
                namesout.write("{}\n".format(p))

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

def rename_nodes_ncbi(tree, verbose=False):
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
    if verbose:
        sys.stderr.write("Traversing the tree to rename the nodes\n")
    for n in tree.traverse("preorder"):
        if n.is_leaf():
            continue
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
                newname = "{} r_{}".format(taxs[w].pop(), w)
                if verbose:
                    sys.stderr.write("Changing name from: {} to {}\n".format(n.name, newname))
                n.name = newname
                break
    return tree

def reroot_tree(tree, verbose=False):
    """
    Reroot the tree between bacteria and archaea.

    This will only work after renaming the leaves on the tree.

    :param tree: the tree
    """
    
    didreroot = False

    if verbose:
        sys.stderr.write("rerooting the tree\n")
    for n in tree.traverse("preorder"):
        childs = n.get_children()
        if verbose:
            cname = ""
            for c in childs:
                cname += "| {} |".format(c.name)
            sys.stderr.write("{}\t{}\t{}\n".format(len(childs), n.name, cname))
        if len(childs) == 2:
            if ("Archaea r_superkingdom" in childs[0].name and "Eukaryota r_superkingdom" in childs[1].name) or ("Archaea r_superkingdom" in childs[1].name and "Eukaryota r_superkingdom" in childs[0].name):
                tree.set_outgroup(n)
                if verbose:
                    sys.stderr.write("Rerooted on {}\n".format(n.name))
                didreroot = True
                break
            if "Bacteria r_superkingdom" in childs[0].name and "Archaea r_superkingdom" in childs[1].name:
                tree.set_outgroup(childs[0])
                if verbose:
                    sys.stderr.write("Rerooted on {}\n".format(childs[0].name))
                didreroot = True
                break
            if "Bacteria r_superkingdom" in childs[1].name and "Archaea r_superkingdom" in childs[0].name:
                tree.set_outgroup(childs[1])
                if verbose:
                    sys.stderr.write("Rerooted on {}\n".format(childs[1].name))
                didreroot = True
                break

    if not didreroot:
        for n in tree.traverse("preorder"):
            if "Bacteria r_superkingdom" in n.name:
                tree.set_outgroup(n)
                if verbose:
                    sys.stderr.write("Rerooted on {} because it is bacteria\n".format(n.name))
                break


    return tree

def write_leaves(tree, outputf):
    """
    Write a list of all the leaves, one line per leaf.
    :param tree: the tree
    :param outputf: the file to write
    :return: 
    """

    with open(outputf, 'w') as out:
        for n in tree.get_leaves():
            out.write("{}\n".format(n.name))
        

def write_tree(tree, outputf):
    """
    Write the tree to a file.
    :param tree: The tree to write
    :param outputf: The output filename
    :return:
    """

    tree.write(outfile=outputf, format=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse a jplacer file')
    parser.add_argument('-j', help='jplacer file', required=True)
    parser.add_argument('-o', help='output file to write the tree to', required=True)
    parser.add_argument('-n', help='names file for nodes that were placed onto the tree. If provided we write out those names')
    parser.add_argument('-l', help='filename for a file with a list of all leaf names of the final tree')
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
    tree = insert_new_nodes(tree, placements, args.n, args.v)
    #tree = rename_leaves_taxids(tree)
    tree = rename_nodes_ncbi(tree, args.v)
    tree = reroot_tree(tree, args.v)
    write_tree(tree, args.o)
    if args.l:
        write_leaves(tree, args.l)
