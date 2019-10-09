"""
Parse a GTO object
"""

import os
import sys
import argparse

from roblib import bcolors
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot a heatmap")
    parser.add_argument('-f', help='gto file', required=True)
    parser.add_argument('-o', help='output file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    t = json.load(open(args.f, 'r'))

    print("{}".format("\n".join(t.keys())))
