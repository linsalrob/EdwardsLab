"""
Generate a random sequence
"""

import os
import sys
import argparse
from random import randint
__author__ = 'Rob Edwards'


def random_sequence(maxlen):
    """
    Generate a random DNA sequence less than maxlen size
    """

    bases = ["A", "G", "C", "T"]
    for i in range(maxlen):
        print(bases[randint(0,3)], end="")
    print()





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-n', type=int, default=1000,
                        help='maximum sequence length (default=1000)')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()
    
    random_sequence(args.n)






