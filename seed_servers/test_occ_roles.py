"""
Test the conversion of roles to pegs
"""

import os
import sys

import argparse
from servers.SAP import SAPserver

def occ_to_roles(roles):
    """
    Convert roles to a dict of roles and the pegs they do

    :param roles:
    :type roles:
    :return:
    :rtype:
    """

    sv = SAPserver()
    result = sv.occ_of_role({'-roles' : roles})
    return result

if __name__ == '__main__':
    roles = ['PTS system, N-acetylglucosamine-specific IIB component (EC 2.7.1.69)', 'Glycerol-3-phosphate dehydrogenase [NAD+] (EC 1.1.1.8)']
    res = occ_to_roles(roles)
    for r in res:
        print(r + "\t" + "\n".join(res[r]))