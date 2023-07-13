"""
Parse a directory of directories of mmseqs eaxy_taxonomies and create a taxonomy table
"""

import os
import sys
import argparse
import gzip

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='mmseqs output directory', required=True)
    parser.add_argument('-e', help='file extension: "_report.')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()
