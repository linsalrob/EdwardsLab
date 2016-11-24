"""
Parse some websites and print some information
"""

import os
import sys
import argparse
from bs4 import BeautifulSoup

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse a directory of websites')
    parser.add_argument('-d', help='directory of files', required=True)
    args = parser.parse_args()

    for f in os.listdir(args.d):
        try:
            aid = f.replace('.aspx', '')
            soup = BeautifulSoup(open(os.path.join(args.d, f), 'r'), 'html.parser')

            print("Organism\t{}".format(soup.title.get_text()))
            print("ID\t{}".format(aid))
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
        except Exception as err:
            sys.stderr.write("There was an error parsing {}\n{}\n. Skipped\n".format(f, err))

