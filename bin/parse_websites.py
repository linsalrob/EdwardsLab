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

            desc = soup.find('meta', attrs={'name':"description"})
            if desc:
                descc = desc['content']
            else:
                descc = "No description"

            if "Designation" in descc:
                parts=re.search('^(.*?)\s+Designation:\s+(.*?)\s+TypeStrain=(\S+)\s+Application:\s*(.*?)$', descc)
                if parts:
                    (org, desig, typestrain, application)=parts.groups()
                    org = org.replace('&reg;', '')
                    org = org.replace('&trade;', '')

                    print("Name\t{}\nDesignation\t{}\nType Strain\t{}\nApplication\t{}\n".format(org, desig, typestrain, application))
                else:
                    sys.stderr.write("Malformed description in {}\n".format(os.path.join(args.d, f)))
            elif descc:
                print("Name\t{}\n".format(descc))
            else:
                sys.stderr.write("No description in {}\n".format(os.path.join(args.d, f)))

            organism = soup.title.get_text()
            organism = organism.replace('\n', '; ')
            organism = organism.replace('\r', '')
            (organism, n) = re.subn(';\s+;', '; ', organism)
            while (n > 0):
                (organism, n) = re.subn(';\s+;', '; ', organism)
            organism = re.sub(';\s+', '; ', organism)
            organism = re.sub('\s+:\s+;\s+', ' : ', organism)
            organism = re.sub('^[\s\;]+', '', organism)
            organism = re.sub('[\s\;]+$', '', organism)



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
                        cell = re.sub('^[\s\;]+', '', cell)
                        cell = re.sub('[\s\;]+$', '', organism)

                    print("{}\t{}".format(header, cell))
            print("//")
        except Exception as err:
            sys.stderr.write("There was an error parsing {}\n{}\n. Skipped\n".format(f, err))

