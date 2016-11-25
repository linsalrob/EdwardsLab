"""
Retrieve just the file sizes for the SRA files
"""

import os
import sys
import argparse
import requests

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Retrieve the file sizes for the SRA data sets')
    parser.add_argument('-f', help='file with tab separated list of SRA IDs and Metagenome Types (e.g. from PARTIE)', required=True)
    args = parser.parse_args()

    with open(args.f, 'r') as f:
        for l in f:
            sraid, type = l.strip().split("\t")
            # According to: https://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.determining_the_location_of
            # The URL is
            # url="ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR[0:3]/srr[0:6]/srr/srr.sra"
            url = "https://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{}/{}/{}/{}.sra".format(sraid[0:3], sraid[0:6], sraid, sraid)
            r = requests.head(url, headers={'Accept-Encoding': 'identity'})
            print("{}\t{}".format(sraid, r.headers['content-length']))



