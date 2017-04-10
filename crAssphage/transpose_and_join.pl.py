"""
Transpose and join a whole lot of files created by coverage_depth.py
"""

import os, sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read a directory of files and join them as a single file")
    parser.add_argument('-d', help='directory of output files from coverage_depth.py', required=True)
    args = parser.parse_args()

    for f in os.listdir(args.d):
        with open(os.path.join(args.d, f), 'r') as fin:
            firstline = True
            for l in fin:
                p=l.strip().split("\t")
                if firstline:
                    firstline = False
                    sys.stdout.write(p[1])
                sys.stdout.write("\t{}".format(p[1]))
            sys.stdout.write("\n")


