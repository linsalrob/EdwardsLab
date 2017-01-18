"""
Export authors from a bib file.

You can use this to create a new conflicts list for 

"""

import os
import sys
import argparse

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    args = parser.parse_args()