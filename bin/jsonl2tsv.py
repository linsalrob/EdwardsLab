"""
Convert a Google JSONL file to tsv.

This format of JSON that comes from big query has one dictionary per line, and is somewhat unique to Google. Note that
the file is not valid JSON format, because each line is an entry.
"""

import os
import sys
import argparse
import json

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert JSON to tab separated format')
    parser.add_argument('-j', '--json', help='JSON input file from a Google bq', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-s', '--separator', help='list separator used to join them (Default = |)',
                        default='|')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    data = []
    allkeys = set()
    with open(args.json, 'r') as jin:
        for l in jin:
            js = json.loads(l)
            thisdata = {}
            for d in js:
                if isinstance(js[d], list):
                    if len(js[d]) == 0:
                        continue
                    if isinstance(js[d][0], dict):
                        # this is a list of key/value pairs!
                        for di in js[d]:
                            mykey = f"{d} : {di['k']}"
                            thisdata[mykey] = di['v']
                    else:
                        thisdata[d] = args.separator.join(js[d])
                else:
                    thisdata[d] = js[d]

            allkeys.update(thisdata.keys())
            data.append(thisdata)

    if 'acc' in allkeys:
        allkeys.discard('acc')
    else:
        print("Could not pop acc to the front of the list")
    if 'jattr' in allkeys:
        allkeys.discard('jattr')

    ak = sorted(list(allkeys))
    ak.insert(0, 'acc')
    ak.append('jattr')
    with open(args.output, 'w') as out:
        # print the header
        print("\t".join(ak), file=out)
        for js in data:
            for k in ak:
                if k not in js:
                    js[k] = ""
                print(f"{js[k]}", end="\t", file=out)
            print(file=out)
