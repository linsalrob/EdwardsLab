"""
Create a color strip file from a tsv file that has the sequence ids.
"""

import os
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create a color strip file")
    parser.add_argument('-f', help='')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()