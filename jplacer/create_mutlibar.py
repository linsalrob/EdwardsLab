"""
Create a multibar file from our labeled file
"""

import os
import sys
import argparse

def read_labels(lf, col, verbose=False):
    """
    Read the labels file and return a dict with tree labels and values
    :param lf: labels file
    :param col: the column to use
    :param verbose: extra output
    :return: a dict
    """

    ret = {}
    with open(lf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if len(p) < col:
                continue
            if not p[col]:
                continue
            ret[p[0]] = p[col]
    return ret




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', help='')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()