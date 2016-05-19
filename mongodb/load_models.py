"""
Create a mongo database if it doesn't exist and load a bunch of data into it.

We need a directory with one or more JSON files in it. We look for JSON on the end of the filename.

e.g. python load_models.py -d /data/Genotype-Phenotype-Modeling/models/Citrobacter/Citrobacter/models/ -n fba_models -c citrobacter
"""

import os
import sys

import argparse
import json
from pymongo import MongoClient

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Load some data from a directory of JSON files")
    parser.add_argument('-d', help='Directory of files', required=True)
    parser.add_argument('-n', help='Database name', required=True)
    parser.add_argument('-c', help='Collection name', required=True)
    args = parser.parse_args()

    client = MongoClient()

    db = client[args.n]
    coll = db[args.c]

    for f in os.listdir(args.d):
        if f.lower().endswith('.json'):
            sys.stderr.write("Loading file " + f + "\n")
            text = json.load(open(os.path.join(args.d, f)))
            obj = {'file_name' : os.path.join(args.d, f), 'content' : text}
            coll.insert(obj)

