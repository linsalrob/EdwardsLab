"""
Parse some websites and print some information
"""

import os
import sys
import re
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

            desc = soup.find('meta', name='description')
            if desc:
                descc = desc['content']
            else:
                descc = "No description"

            parts=re.search('^(.*?)\s+Designation:\s+(.*?)\s+TypeStrain=(\S+)\s+Application:\s+(.*?)$', descc)
            (org, desig, typestrain, application)=parts.groups()
            print("Name\t{}\nDesignation\t{}\nType Strain\t{}\nApplication\t{}\n".format(org, desig, typestrain, application))

            organism = soup.title.get_text()
            organism = organism.replace('\n', '; ')
            organism = organism.replace('\r', '')
            (organism, n) = re.subn(';\s+;', '; ', organism)
            while (n > 0):
                (organism, n) = re.subn(';\s+;', '; ', organism)
            organism = re.sub(';\s+', '; ', organism)
            organism = re.sub('\s+:\s+;\s+', ' : ', organism)
            organism = organism.strip()



            print("Organism\t{}".format(organism))

            print("ID\t{}".format(aid))
            for table in soup.find_all('table', class_="fulllist"):
                for row in table.find_all('tr'):
                    header = row.find('th').get_text()
                    cell = row.find('td').get_text()

                    if header:
                        header = header.strip()
                        # header = header.encode('ascii', 'ignore')

                    if cell:
                        # cell = cell.encode('ascii', 'ignore')
                        cell = cell.replace('\n', '; ')
                        cell = cell.replace('\r', '')
                        (cell, n) = re.subn(';\s+;', '; ', cell)
                        while (n > 0):
                            (cell, n) = re.subn(';\s+;', '; ', cell)
                        cell = re.sub(';\s+', '; ', cell)
                        cell = re.sub('\s+:\s+;\s+', ' : ', cell)
                        cell = cell.strip()

                    print("{}\t{}".format(header, cell))
            print("//")
        except Exception as err:
            sys.stderr.write("There was an error parsing {}\n{}\n. Skipped\n".format(f, err))

