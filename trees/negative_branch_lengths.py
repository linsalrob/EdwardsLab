"""
Accoring to http://www.icp.ucl.ac.be/~opperd/private/neighbor.html:

As the neighbor-joining algorithm seeks to represent the data in the form of an additive tree, it can assign a negative
length to the branch. Here the interpretation of branch lengths as an estimated number of substititions gets into
difficulties. When this occurs it is adviced to set the branch length to zero and transfer the difference to the
adjacent branch length so that the total distance between an adjacent pair of terminal nodes remains unaffected.
This does not alter the overall topology of the tree (Kuhner and Felsenstein, 1994).

We will identify negative branch lengths, set them to 0 and transfer the remainder.
"""

import os
import sys
from roblib import Newick_Tree
import argparse


def find_negative(node):
    """
    Find a negative distance
    :param node: the node we're at
    :return:
    """
    if node.distance < 0:
        sys.stderr.write("Negative distance: {}\n".format(node.distance))

    if node.left:
        find_negative(node.left)
    if node.right:
        find_negative(node.right)



def correct_negative(node):
    """
    Correct a negative distance at the left or right children of node
    :param node: the node we're at
    :return:
    """

    if node.left and node.left.distance < 0:
        if node.right:
            node.right.distance -= node.left.distance  # note that this adds the abs value of dist to node.right
        else:
            node.distance -= node.left.distance
        node.left.distance = 0

    if node.right and node.right.distance < 0:
        if node.left:
            node.left.distance -= node.right.distance
        else:
            node.distance -= node.right.distance
        node.right.distance = 0

    if node.left:
        correct_negative(node.left)
    if node.right:
        correct_negative(node.right)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse a tree and correct negative branch lengths')
    parser.add_argument('-t', help='tree file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    tre = []
    with open(args.t, 'r') as f:
        for l in f:
            tre.append(l.strip())

    root = Newick_Tree().parse(''.join(tre))
    if args.v:
        sys.stderr.write("NEGATIVE NODES:\n")
        find_negative(root)
    correct_negative(root)
    Newick_Tree().print_tree(root)


