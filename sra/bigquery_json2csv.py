"""
Parse the output from a bigquery output and convert to CSV
"""
import json
import os
import sys
import argparse
import json

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-b', help='big query output file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    keys = set()
    allattr = set()
    acc = {}
    attr = {}
    with open(args.b, 'r') as f:
        for l in f:
            l = l.strip()
            j = json.loads(l)
            acc[j['acc']] = {}
            attr[j['acc']] = {}
            for k in j:
                if 'jattr' == k:
                    continue
                if 'attributes' == k:
                    for piece in j['attributes']:
                        allattr.add(piece['k'])
                        attr[j['acc']][piece['k']] = piece['v']
                    continue
                acc[j['acc']][k] = j[k]
                keys.add(k)

    allkeys = sorted(keys)
    attrs = sorted(allattr.difference(keys))
    with open(args.o, 'w') as out:
        print("\t".join(["Accession"] + allkeys + attrs), file=out)
        for a in sorted(acc.keys()):
            print(a, file=out, end="")
            for k in allkeys:
                if k in acc[a]:
                    print(f"\t{acc[a][k]}", file=out, end="")
                else:
                    print("\t", file=out, end="")
            for m in attrs:
                if m in attr[a]:
                    print(f"\t{attr[a][m]}", file=out, end="")
                else:
                    print("\t", file=out, end="")

            print(file=out)


