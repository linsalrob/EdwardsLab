"""
A very simple JSON validator. We read the JSON file in and print it out using pprint.
"""

import sys
import json 
import pprint
__author__ = 'Rob Edwards'


h = "A very simple JSON validator. We read the JSON file in and print it out using pprint.\n"
h += f"\nUsage: {sys.argv[0]} <json file>\n"


if len(sys.argv) < 2:
    sys.exit(h)


pp = pprint.PrettyPrinter(indent=4)

with open(sys.argv[1], 'r') as f:
    j = json.load(f)

pp.pprint(j)

