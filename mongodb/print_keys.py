"""
Print the keys from a file
"""

import os
import sys

import argparse
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print all the keys in a file")
    parser.add_argument('-f', help='Files to print keys from', required=True)
    args = parser.parse_args()

    data = json.load(open(args.f, 'r'))
    for k in data:
        print(k)