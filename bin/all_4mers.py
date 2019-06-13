"""
print all 4 mers
"""

import os
import sys
import argparse

from roblib import bcolors
import itertools

alphabet = ['A', 'C', "G", 'T']
c=0
for a in itertools.product(alphabet, repeat=4):
    if (c == 16):
        print()
        c=0
    print("".join(a), end=' ')
    c+=1
print()