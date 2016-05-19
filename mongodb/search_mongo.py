"""
Search a mongo database loaded with load_models.py
'file_name' : '/data/Genotype-Phenotype-Modeling/models/Citrobacter/Citrobacter/models/C.sedlakii_gf_draft_ArgonneLB.json'

"""


import argparse

from pymongo import MongoClient

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Search a mongo database")
    parser.add_argument('-n', help='Database name', required=True)
    parser.add_argument('-c', help='Collection name', required=True)
    parser.add_argument('-k', help='Key to find', required=True)
    parser.add_argument('-v', help='Value for key', required=True)
    args = parser.parse_args()

    client = MongoClient()

    db = client[args.n]
    coll = db[args.c]

    for cursor in coll.find({args.k : args.v}):
        # this is a cursor to the document
        print(cursor['file_name'])
