import os
import sys
import argparse


def genera(envF):
    """
    Read genera per environment
    :param envF:
    :return:
    """

    with open(envF, 'r') as f:







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate shared genera per env")
    parser.add_argument('-f', help='file of environments and genera')
    args = parser.parse_args()