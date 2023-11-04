#!/usr/bin/env python

"""
Read a tsv file and transpose it. 
"""

import os
import sys
import argparse
import pandas as pd
from roblib import bcolors
__author__ = 'Rob Edwards'





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-s', '--sep', help='separator (default=tab)', default="\t")
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.v:
        print(f"{bcolors.GREEN}Reading {args.input}{bcolors.ENDC}")
    df = pd.read_csv(args.input, sep=args.sep)
    dft = df.T
    if args.v:
        print(f"{bcolors.GREEN}Writing {args.output}{bcolors.ENDC}")
    dft.to_csv(args.output, sep=args.sep)
    



