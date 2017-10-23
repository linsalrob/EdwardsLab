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

We use an idea inspired by the ete3 team: https://gist.github.com/jhcepas/279f9009f46bf675e3a890c19191158b :

For each node find its path to the root.

e.g.

A -> A, y, z
E -> E, w, x,z

and make these orderless sets. Then we XOR the two sets to only find the elements
that are in one or other sets but not both. In this case A, E, y, x, w.

The distance between the two nodes is the sum of the distances from each of those nodes
to the parent

One more optimization: since the distances are symmetric, and distance to itself is zero
we user itertools.combinations rather than itertools.permutations. This cuts our computes from theta(n^2)
1/2n^2 - n (= O(n^2), which is still not great, but in reality speeds things up for large trees).


"""

import os
import sys
import argparse
from itertools import combinations
from ete3 import Tree


def make_matrix(treefile):
    """
    Create a matrix from a tree file
    :param treefile:
    :return:
    """

    tree = Tree(treefile)

    leaves = tree.get_leaves()
    paths = {x:set() for x in leaves}

    # get the paths going up the tree
    # we get all the nodes up to the last one and store them in a set
    sys.stderr.write("Precalculating distances\n")
    for n in leaves:
        if n.is_root():
            continue
        movingnode = n
        while not movingnode.is_root():
            paths[n].add(movingnode)
            movingnode = movingnode.up

    # now we want to get all pairs of nodes using itertools combinations. We need AB AC etc but don't need BA CA

    leaf_distances = {x.name:{} for x in leaves}

    sys.stderr.write("Iterating over the leaves\n")
    for (leaf1, leaf2) in combinations(leaves, 2):
        # figure out the unique nodes in the path
        uniquenodes = paths[leaf1] ^ paths[leaf2]
        distance = sum(x.dist for x in uniquenodes)
        leaf_distances[leaf1.name][leaf2.name] = leaf_distances[leaf2.name][leaf1.name] = distance

    allleaves = sorted(leaf_distances.keys())
    sys.stdout.write("\t".join([""] + allleaves) + "\n")
    for n in allleaves:
        sys.stdout.write(n + "\t")
        for m in allleaves:
            if m == n:
                sys.stdout.write("0\t")
            else:
                sys.stdout.write("{}\t".format(leaf_distances[n][m]))
        sys.stdout.write("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a tree into a distance matrix')
    parser.add_argument('-t', help='Tree file', required=True)
    args = parser.parse_args()

    make_matrix(args.t)