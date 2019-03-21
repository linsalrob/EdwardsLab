"""
Parse a directory with html files with all the student theses in them
"""

import os
import argparse
from bs4 import BeautifulSoup
import re
import sys

def parse_thesis(f, verbose=False):
    """
    Parse a single thesis and return a dict of the information
    :param f: The file name to parse
    :param verbose: more output
    :return: a tuple of the RedID and a dict of information
    """

    if verbose:
        sys.stderr.write(f"Parsing {f}\n")

    soup = BeautifulSoup(open(f, 'r'), 'html.parser')
    data = {}

    # get the student information
    for s in soup.body.find_all(text=re.compile('RedID')):
        t = s.find_parent('table')
        for row in t.find_all("tr"):
            r = [re.sub('\s+', ' ', cell.get_text(strip=True)) for cell in row.find_all("th")]
            s = [re.sub('\s+', ' ', cell.get_text(strip=True)) for cell in row.find_all("td")]
            data.update(dict(zip(r, s)))

    if verbose:
        sys.stderr.write(f"\t{data['First Name']} {data['Last Name']} ")

    # thesis status
    for s in soup.body.find_all(text=re.compile('Submitted')):
        t = s.find_parent('table')
        rows = t.find_all("tr")
        r = [re.sub('\s+', ' ', cell.get_text(strip=True)) for cell in rows[0].find_all("th")]
        s = [re.sub('\s+', ' ', cell.get_text(strip=True)) for cell in rows[1].find_all("td")]
        l = [re.sub('\s+', ' ', cell.get_text(strip=True)) for cell in rows[2].find_all("td")][1]
        data.update({"Thesis Title": l})
        data.update(dict(zip(r, s)))

    # thesis committee
    committee = {}
    for s in soup.body.find_all(text=re.compile('Chair')):
        t = s.find_parent('table')
        for row in t.find_all("tr"):
            cells = row.find_all("td")
            what = cells[0].get_text(strip=True)
            who = cells[1].get_text(strip=True)
            committee.update(dict(zip([what], [who])))
    data.update(committee)

    if verbose:
        sys.stderr.write(f" Chair: {data['Chair']}\n")

    return data['RedID'], data

def print_data(data, outputfile, verbose=False):
    """
    Print the data to a file
    :param data: The dict of information, with red id as a key to a dict of all information
    :param outputfile: The file to write
    :param verbose: more output
    :return: nothing
    """

    cols = set()

    for r in data:
        cols.update(data[r].keys())

    wanted = ['RedID', 'First Name', 'Last Name', 'Primary Major', 'Primary Major Name', 'Degree Objective', 'Degree Conferred', 'Approved', 'Approved date', 'Submitted', 'Published', 'Chair', '2nd Member', '3rd Member', '4th Member', '5th Member']
    others = list(cols.difference(wanted))


    with open(outputfile, 'w') as out:
        out.write("redid\t" + "\t".join(wanted + others) + "\n")
        for redid in data:
            out.write(f"{redid}")
            for w in wanted:
                out.write("\t" + data[redid].get(w, ""))
            for w in others:
                out.write("\t" + data[redid].get(w, ""))
            out.write("\n")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse a directory of theses html pages")
    parser.add_argument('-d', help='Directory of theses html pages to parse', required=True)
    parser.add_argument('-o', help='Output file to write', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    students = {}
    for f in os.listdir(args.d):
        r, d = parse_thesis(os.path.join(args.d, f), args.v)
        students[r] = d

    print_data(students, args.o, args.v)