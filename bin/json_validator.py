"""

"""

import os
import sys
import json 
import pprint
import argparse
from roblib import bcolors
__author__ = 'Rob Edwards'

pp = pprint.PrettyPrinter(indent=4)
f='process_metagenomes.json'

with open(f, 'r') as fin:
    j = json.load(fin)

pp.pprint(j)

