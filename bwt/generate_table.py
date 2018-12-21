"""
Generate a table with all the BWT rotations
"""

import os
import sys
import argparse

__author__ = 'Rob Edwards'

def printarr(arr):
    """
    print the array
    :param arr:
    :return:
    """
    print("\n".join(arr))


def rotate(s):
    """
    Rotate the string
    :param s:
    :return:
    """

    d=[]
    d.append(s)
    for i in reversed(range(1, len(s))):
        d.append(s[i:] + s[:i])
    return d

def count_ends(d):
    """
    count number of consecutive letters
    :param d:
    :return:
    """
    con=0
    for i in range(len(d)-1):
        if d[i][-1] == d[i+1][-1]:
            con+=1
    print("{} consecutive letters".format(con))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate all the rotations of a sequence')
    parser.add_argument('-s', help='string to rotate', required=True)
    parser.add_argument('-b', help="character to start the string with (defualt = no character)", default=None)
    parser.add_argument('-e', help='Character to end the string with (default = $)', default='$')
    parser.add_argument('-c', help='change case', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    s = args.s + args.e
    if args.b:
        s = args.b + s

    if args.c:
        s=s.lower()

    d = rotate(s)
    printarr(d)
    count_ends(d)

    print("\nAfter sorting:")
    d = sorted(d)
    printarr(d)

    count_ends(d)

