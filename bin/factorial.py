"""
Calculate some factorials
"""

import os
import sys


def factorial(n):
    if n == 2:
        return 2
    return n * factorial(n-1)


for n in range(2,100):
    if factorial(n) > 100000 ** n:
        print("{}\tFACTORIAL".format(n))
    else:
        print("{}\texponent".format(n))
#    print("{}\t{}\t{}".format(n, factorial(n), 100000 ** n))