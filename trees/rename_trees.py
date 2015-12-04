import os
import sys
import argparse
import os
import string
import sys

__author__ = 'Rob Edwards'

"""
Rename tree leaves and print out the tree
"""

class Node(object):
    """A node object"""

    def __init__(self, id):
        self.id = id
        self.left = None
        self.right = None
        self.parent = None
        self.distance = ""
        self.name = ""
        self.side = None


class Tree(object):
    def __init__(self):
        pass

    def count_nodes(self, root):
        """ Count the number of nodes in the tree
        :param root: The root node
        :type root: Node
        :return: The number of nodes
        :rtype: int
        """

        def count(node):
            c = 1
            if node.left:
                c += count(node.left)
            if node.right:
                c += count(node.right)
            return c

        return count(root)

    def parse(self, tree):
        def process_tree(treestr, pos, node):
            # sys.stderr.write("At pos {} tree has depth {}\n".format(pos, Tree().count_nodes(root)))

            if treestr[pos] == '(':
                pos += 1
                # sys.stderr.write("LEFT: When adding node {} tree has depth {}\n".format(pos, Tree().count_nodes(root)))
                newnode = Node(pos)
                newnode.parent = node
                node.left = newnode
                newnode.side = "Left"
                # sys.stderr.write("ADDED NODE {} to the LEFT\n".format(newnode.id))
                return process_tree(treestr, pos, newnode)
            elif treestr[pos] == ',':
                pos += 1
                # sys.stderr.write("RIGHT: When adding node {} tree has depth {}\n".format(pos, Tree().count_nodes(root)))
                newnode = Node(pos)
                if node.parent.right:
                    newnode = node.parent
                else:
                    newnode.parent = node.parent
                    node.parent.right = newnode
                    newnode.side = 'Right'
                # sys.stderr.write("ADDED NODE {} to the RIGHT\n".format(newnode.id))
                return process_tree(treestr, pos, newnode)
            elif treestr[pos] == ')':
                pos += 1
                if pos >= len(treestr):
                    return
                while treestr[pos] in string.ascii_letters or treestr[pos] == '_' or treestr[pos] in string.digits or \
                                treestr[pos] in '-':
                    node.name += treestr[pos]
                    pos += 1
                return process_tree(treestr, pos, node.parent)
            elif treestr[pos] == ':':
                pos += 1
                try:
                    while treestr[pos] in string.digits or treestr[pos] == '.' or treestr[pos] in '-':
                        node.distance += treestr[pos]
                        pos += 1
                except TypeError:
                    raise TypeError(
                        "TypeError: CANNOT ADD {} at POS {} to {} in node {}\n".format(treestr[pos], pos, node.distance,
                                                                                       node.id))
                # sys.stderr.write("Set NODE {} dist to {}\n".format(node.id, node.distance))
                node.distance = float(node.distance)
                return process_tree(treestr, pos, node)
            else:
                while treestr[pos] in string.ascii_letters or treestr[pos] == '_' or treestr[pos] in string.digits or \
                                treestr[pos] in '-':
                    node.name += treestr[pos]
                    pos += 1
                # sys.stderr.write("When adding node {} tree has depth {}\n".format(node.name, Tree().count_nodes(root)))
                return process_tree(treestr, pos, node)

        parent = Node("root")
        root = parent
        pos = 0
        tree = tree.rstrip(';')
        pos = process_tree(tree, pos, parent)
        return parent

    def print_tree(self, root):
        """
        Print out a tree in newick format.

        :param root: The root node of the tree
        :type root: Node
        :return:
        :rtype:
        """

        def process_child(node):
            toreturn = ''
            if node.left or node.right:
                toreturn = '('
            if node.left and node.right:
                toreturn += process_child(node.left) + "," + process_child(node.right)
            elif node.left:
                toreturn += process_child(node.left)
            elif node.right:
                toreturn += process_child(node.right)
            if node.left and node.right and node.name:
                # the root node??
                toreturn += ','
            elif node.left or node.right:
                toreturn += ")"
            if node.name:
                toreturn += node.name
            toreturn += ":{}".format(node.distance)
            return toreturn

        print(process_child(root) + ");")

    def clean_name(self, name):
        """
        Just clean out non-allowable characters in the name

        :param name: the new name
        :type name: str
        :return: the revised name
        :rtype: str
        """

        allowable = set(string.ascii_letters)
        allowable.update(set(string.digits))
        allowable.update({'_','-',':'})
        name.replace(' ', '_')
        return filter(lambda x: x in allowable, name)

    def rename_nodes(self, node, idmap):
        """
        Rename the nodes of a tree based on id map

        :param root: the root node of the tree
        :type root: Node
        :param idmap: the id map
        :type idmap: dict
        :return: the new root node
        :rtype: Node
        """

        if node.name and node.name in idmap:
            node.name = Tree().clean_name(idmap[node.name])
        if node.left:
            node.left = Tree().rename_nodes(node.left, idmap)
        if node.right:
            node.right = Tree().rename_nodes(node.right, idmap)

        return node


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse a tree')
    parser.add_argument('-t', help='tree file', required=True)
    parser.add_argument('-i', help='id map file', required=True)
    args = parser.parse_args()

    idmap = {}
    with open(args.i, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            idmap[p[0]]=p[1].split()[0]


    tre = []
    with open(args.t, 'r') as f:
        for l in f:
            tre.append(l.strip())

    root = Tree().parse(''.join(tre))
    root = Tree().rename_nodes(root, idmap)
    Tree().print_tree(root)

