"""
Search for something in the database

"""

import sys
import os

import json
from pymongo import MongoClient

client = MongoClient()

db = client['fba_models']
coll = db['citrobacter']

for cursor in coll.find({ "media_ref": "5067/41/1" }):
    # examples from the models data;
    """
    "content.id" :  "contig00112.fbamdl3.fba.90"

    "content.objectiveValue"  : { '$gt' : 10 }
    "media_ref":"5067/41/1"

    """
    # this is a cursor to the document
    print(cursor['file_name'])


