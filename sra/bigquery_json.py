"""
Parse the output from a bigquery output from JSON and save some stuff
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
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    with open(args.b, 'r') as f:
        for l in f:
            l = l.strip()
            j = json.loads(l)
            ll = None
            if 'attributes' in j:
                for a in j['attributes']:
                    if a['k'] == 'lat_lon_sam':
                        ll = a['v']
            gls = None
            if "geo_loc_name_sam" in j:
                if len(j['geo_loc_name_sam']) > 0:
                    gls = j['geo_loc_name_sam'][0]
            print(f"{j['acc']}\t{gls}\t{ll}")
