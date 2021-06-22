"""
Read some json files and print all keys
"""

import os
import sys
import argparse

from roblib import bcolors
import json


__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='json file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    akeys = set()
    t = json.load(open(args.f, 'r'))
    for k in t:
        akeys.update(k.keys())

    print("{}".format("\n".join(akeys)))
