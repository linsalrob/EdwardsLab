"""
Start with a tree file, use ete3 to create a pairwise distance for all nodes. Basically this is the
distance matrix but as tuples.

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


def make_dists(treefile, printone, verbose):
    """
    Create pairwise distances from a tree file
    :param treefile: the tree file to parse
    :param printone: if true we only print one copy of the pair (ie. A -> B). If false we print A->B and B->A
    :param verbose: make some additional output
    :return:
    """

    tree = Tree(treefile)

    leaves = tree.get_leaves()
    paths = {x:set() for x in leaves}

    # get the paths going up the tree
    # we get all the nodes up to the last one and store them in a set
    if verbose:
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

    if verbose:
        sys.stderr.write("Iterating over the leaves\n")
    for (leaf1, leaf2) in combinations(leaves, 2):
        # figure out the unique nodes in the path
        uniquenodes = paths[leaf1] ^ paths[leaf2]
        distance = sum(x.dist for x in uniquenodes)
        if printone:
            if leaf1.name < leaf2.name:
                print("{}\t{}\t{}".format(leaf1.name, leaf2.name, distance))
            else:
                print("{}\t{}\t{}".format(leaf2.name, leaf1.name, distance))
        else:
            print("{}\t{}\t{}".format(leaf1.name, leaf2.name, distance))
            print("{}\t{}\t{}".format(leaf2.name, leaf1.name, distance))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a tree into a distance matrix')
    parser.add_argument('-t', help='Tree file', required=True)
    parser.add_argument('-p', help='Print one direction (A->B). Default is to print A->B and B->A', action='store_true')
    parser.add_argument('-v', help='Verbose output. (Mostly progress)', action='store_true')
    args = parser.parse_args()

    make_dists(args.t, args.p, args.v)