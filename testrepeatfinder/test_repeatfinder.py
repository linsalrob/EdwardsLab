"""
Test the implementation of the repeatfinder extension
"""

import os
import sys
import argparse

from roblib import bcolors

import RobRepeatFinder
import pprint

s = "TTTTTTTTTTTTagcaTTTTTTTTTTTT"
print(f"s: {s}")
r = RobRepeatFinder.repeatFinder(s, 0)
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(r)

