"""
Create a mongo database if it doesn't exist and load a bunch of data into it.

We need a directory with one or more JSON files in it. We look for JSON on the end of the filename.

e.g. python load_models.py -d /data/Genotype-Phenotype-Modeling/models/Citrobacter/Citrobacter/models/ -n fba_models -c citrobacter
"""


import sys
import os

import json
from pymongo import MongoClient

client = MongoClient()

db = client['fba_models']
coll = db['citrobacter']


for f in os.listdir('/data/Genotype-Phenotype-Modeling/models/Citrobacter/Citrobacter/models/'):
    if f.lower().endswith('.json'):
        sys.stderr.write("Loading file " + f + "\n")
        text = json.load(open(os.path.join(args.d, f)))
        obj = {'file_name': os.path.join(args.d, f), 'content': text}
        coll.insert(obj)

