"""
PArse the dates and count number per month
"""

import os
import sys
import argparse
import datetime

from roblib import bcolors

def get_classification(srafile="/home/redwards/GitHubs/partie/SRA_Metagenome_Types.tsv"):

    cl = {}
    with open(srafile, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            cl[p[0]] = p[1]
    return cl


def parse_dates(datefile="/home/redwards/GitHubs/partie/SRA_Metagenome_Sizes.tsv"):
    """
    Parse and count dates. We use the ReleaseDate information
    :param datefile:
    :return:
    """

    classification = get_classification()

    count = {}
    with open(datefile, 'r') as f:
        header = f.readline()
        p = header.rstrip().split("\t")
        rdcol = None
        if ("ReleaseDate" in p):
            rdcol = p.index("ReleaseDate")
        else:
            sys.stderr.write(f"{bcolors.RED}FATAL:{bcolors.ENDC} {bcolors.BLUE}ReleaseDate{bcolors.ENDC} was not a header in the first row of {datefile}")
            sys.exit(-1)

        for l in f:
            p = l.rstrip().split("\t")
            if (p[0] not in classification):
                sys.stderr.write(f"{p[0]}\n")
                continue
            if (classification[p[0]] != "WGS"):
                continue
            if len(p) < rdcol:
                continue
            date_time_obj = datetime.datetime.strptime(p[rdcol], '%Y-%m-%d %H:%M:%S')
            ym = f"{date_time_obj.month}/{date_time_obj.year}"
            count[ym] = count.get(ym, 0) + 1

    for ym in count:
        print(f"{ym}\t{count[ym]}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', help='file with dates to parse')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    parse_dates()