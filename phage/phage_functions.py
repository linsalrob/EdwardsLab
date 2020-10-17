"""

Read a file with a list of functions (one per line) and count how many are

phage
not phage
hypothetical

"""

import os
import sys
import argparse
from PhiSpyModules import is_phage_func, is_unknown_func
from roblib import is_hypothetical

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def run(filename, verbose=False):

    hypo = 0
    phage = 0
    notp  = 0

    with open(filename, 'r') as f:
        for l in f:
            fn = l.strip()
            if is_unknown_func(fn) or is_hypothetical(fn):
                hypo += 1
            elif is_phage_func(fn):
                phage += 1
            else:
                notp += 1

    total = hypo + phage + notp

    print(f"Phage proteins: {phage} ({(phage/total)*100:2} %)")
    print(f"Hypothetical proteins: {hypo} ({(hypo/total)*100:2} %)")
    print(f"Not phage proteins {notp} ({(notp/total)*100:2} %)")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='filename', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    run(args.f, args.v)
