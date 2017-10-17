"""
Start with a tree file, use ete3 to create a cophenetic distance matrix

if we have a tree like

                           ----A
              _____________|y
              |            |
              |            ----B
      ________|z
              |            ----C
              |            |
              |____________|x     -----D
                           |      |
                           |______|w
                                  |
                                  |
                                  -----E

Where w,x,y,z are internal nodes.
d(A,B) = d(y,A) + d(y,B)
and
d(A, E) = d(z,A) + d(z, E) = {d(z,y) + d(y,A)} + {d(z,x) + d(x,w) + d(w,E)}

The approach that we use is to start at each leave and calculate the distance from that leave to all the parental
internal nodes. Then we start at the first leave and iterate over its parents, getting all the leaf nodes of that parent.
If we don't have a distance we add the new distance. If we do, we skip it. We need to do it this way, otherwise we
could add the incorrect distance. e.g. when we start at D we move to parent (w) and add d(D,E) = d(D,E). Then we move to parent (x)
and add d(D,C) = d(C,D). At this point we don't readd d(D,E) via w and x as it would be wrong.


"""

import os
import sys
import argparse
from ete3 import Tree

def enumerate_leafs(tree, leaves, distances, leaf_distances):
    sys.stderr.write("Enumerating\n")
    for i, node in enumerate(leaves):
        j = i + 1
        sys.stderr.write("{}\n".format(i))
        while j < len(leaves):
            common = tree.get_common_ancestor(node, leaves[j])
            if common not in distances[node]:
                sys.exit("The common node betwee {} and {} ({}) is not in distances".format(node.name, leaves[j].name,
                                                                                            common.name))
            if common not in distances[leaves[j]]:
                sys.exit("The common node betwee {} and {} ({}) is not in distances".format(node.name, leaves[j].name,
                                                                                            common.name))

            leaf_distances[node.name][leaves[j].name] = distances[node][common] + distances[leaves[j]][common]
            leaf_distances[leaves[j].name][node.name] = distances[node][common] + distances[leaves[j]][common]
            j += 1
    return leaf_distances


def make_matrix(treefile):
    """
    Create a matrix from a tree file
    :param treefile:
    :return:
    """

    tree = Tree(treefile)

    distances = {}
    leaves = []
    internalcounter = 0
    for node in tree.traverse():
        if node.is_leaf():
            leaves.append(node)

    # precalculate all the distances up the tree
    sys.stderr.write("Precalculating distances\n")
    for n in leaves:
        if n.is_root():
            continue
        distances[n]={}
        movingnode = n
        dist = 0
        while not movingnode.is_root():
            dist += movingnode.dist
            movingnode = movingnode.up
            distances[n][movingnode] = dist


    sys.stderr.write("\tDone\n")
    leaf_distances = {}
    for l in leaves:
        leaf_distances[l.name] = {}

    # slow way:
    # enumerate_leafs(tree, leaves, distances, leaf_distances)

    # note should be able to do this by visiting a node and going up the tree
    # ignoring any pairs of nodes that we have seen before.
    # should be quicker than getting lca
    sys.stderr.write("Iterating over the nodes\n")
    for node in leaves:
        parent = node
        while not parent.is_root():
            parent = parent.up
            for leaf in parent:
                if leaf.name in leaf_distances[node.name]:
                    continue
                if not leaf.is_leaf():
                    continue

                if leaf.name == node.name:
                    leaf_distances[node.name][leaf.name] = leaf_distances[leaf.name][node.name] = 0
                    continue
                leaf_distances[node.name][leaf.name] = leaf_distances[leaf.name][node.name] = distances[node][parent] + distances[leaf][parent]





    sys.stdout.write("\t" + "\t".join([n.name for n in leaves]) + "\n")
    for n in leaves:
        sys.stdout.write(n.name + "\t")
        for m in leaves:
            sys.stdout.write("{}\t".format(leaf_distances[n.name][m.name]))
        sys.stdout.write("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a tree into a distance matrix')
    parser.add_argument('-v', help='verbose output', action='store_true')
    parser.add_argument('-t', help='Tree file', required=True)
    args = parser.parse_args()

    make_matrix(args.t)