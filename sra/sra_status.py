"""
Report the status of an SRA run
"""

import os
import sys
import argparse
import requests
import json


def get_status(runid, verbose, url='https://www.ncbi.nlm.nih.gov/Traces/sra/status/srastatrep.fcgi/acc-mirroring?acc='):
    """
    Get the status of the run
    :param runid: run id to get
    :param verbose: more output
    :param url: the base url to append the status to
    :return:
    """

    req = url + runid
    if args.v:
        sys.stderr.write("Getting {}\n".format(req))
    r = requests.get(req)
    if args.v:
        sys.stderr.write("Status: {}\n".format(r.status_code))
    d = json.loads(r.text)
    return d

def print_full(runid, data, verbose):
    """
    Print the full output
    :param runid: the run id to check
    :param data: the json object
    :param verbose: more output
    :return:
    """

    for r in data['rows']:
        if r[0] != runid:
            sys.stderr.write("Expected an accession of {} but found {}\n".format(runid))
        for i, j in zip(data['column_names'], r):
            print("{}\t{}".format(i, j))

def print_status(runid, data, verbose):
    """
    Print the status of the run
    :param runid: the run id to check
    :param data:
    :param verbose:
    :return:
    """

    s = data['column_names'].index('Status')
    for r in data['rows']:
        if r[0] != runid:
            sys.stderr.write("Expected an accession of {} but found {}\n".format(runid))
        print("{}\t{}".format(r[0], r[s]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Return the status of an SRA run")
    parser.add_argument('-r', help='SRA Run ID', required=True)
    parser.add_argument('-f', help="full output", action='store_true')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    data = get_status(args.r, args.v)
    if args.f:
        print_full(args.r, data, args.v)
    else:
        print_status(args.r, data, args.v)