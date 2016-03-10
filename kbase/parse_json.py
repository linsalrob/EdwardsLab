"""
Parse the JSON file and print some stuff out
"""

import os
import sys

import argparse
import json

import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse the JSON file downloaded from KBase")
    parser.add_argument('-f', help='JSON file', required=True)
    args = parser.parse_args()

    data = json.load(open(args.f))

    # print all the keys
    print("\n".join(data.keys()))

    # print keys associated with features
    feats = data['features']

    print("There are " + str(len(data['contig_lengths'])) + " contigs")

    sys.exit(0)

    for f in feats:
        #print("\t".join([f['id'], f['function']]))
        for p in f['function'].split(' ; '):
            m = re.match('\s*[\d\-\.]+$', p)
            if m and m.end() == len(p):
                print("\t".join([f['id'], 'EC ' + p]))
           # else:
           #     print("\t".join([f['id'], p]))


