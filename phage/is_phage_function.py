"""
Note that this uses code from phispy!!

Test whether functions are phages .... or not!
"""

import os
import sys
import argparse

from PhiSpyModules import is_phage_func, is_unknown_func

def is_phage_hypo(f):
    with open(f, 'r') as fin:
        for l in fin:
            p = l.strip().split("\t")
            if is_phage_func(p[0]):
                p.append(1)
            else:
                p.append(0)
            if is_unknown_func(p[0]):
                p.append(1)
            else:
                p.append(0)
            print("\t".join(map(str, p)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='file', required=True)
    args = parser.parse_args()

    is_phage_hypo(args.f)
