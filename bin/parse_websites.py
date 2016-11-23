"""
Parse some websites and print some information
"""

import os
import sys
import argparse
from bs4 import BeautifulSoup

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description=' ')
    # parser.add_argument('-f', help='input file', required=True)
    # args = parser.parse_args()

    # replace file with args.f

    file = "/home/redwards/Desktop/ATCC/10004.aspx"

    soup = BeautifulSoup(open(file, 'r'), 'html.parser')

    print("Organism\t{}".format(soup.title.get_text()))
    for table in soup.find_all('table', class_="fulllist"):
        for row in table.find_all('tr'):
            header = row.find('th').get_text()
            cell = row.find('td').get_text()

            if header:
                header = header.strip()
                # header = header.encode('ascii', 'ignore')

            if cell:
                cell = cell.strip()
                # cell = cell.encode('ascii', 'ignore')

            print("{}\t{}".format(header, cell))
    print("//")