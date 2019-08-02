"""
Read some json files and print all keys
"""

import os
import sys
import argparse
import json


__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='json file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    akeys = set()
    t = json.load(open(args.f, 'r'))

    print("{}".format("\n".join(t.keys())))
