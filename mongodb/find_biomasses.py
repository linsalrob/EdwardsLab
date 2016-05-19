"""
Find biomasses greater than a value
"""

import os
import sys

from pymongo import MongoClient

databasename = 'fba_models'
collectionname = 'citrobacter'

value = 40


client = MongoClient()

coll = client[databasename][collectionname]

for cursor in coll.find({"content.biomasses.biomasscompounds.coefficient": { '$gt': 35 }}):
    print(cursor['file_name'])