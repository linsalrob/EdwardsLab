"""
Convert an xml file to a csv file. We're looking for specific things here
"""

import os
import sys
import argparse
from bs4 import BeautifulSoup
import csv
from dateutil import parser


def read_xml(f):
    """
    Read the xml file
    :param f: File to read
    :return:
    """
    soup = BeautifulSoup(open(f, 'r'), 'xml')

    data  = {}
    for d in soup.find_all("Device"):
        data[d['name']] = []
        for p in d.find_all('Position'):
            g = p.find('GpsAt')
            l = p.find('Latitude')
            n = p.find('Longitude')
            t = parser.parse(g.text)

            timestamp = "{:0>2d}{:0>2d}{:0>2d}{:0>2d}{:0>2d}".format(t.year - 2000, t.month, t.day, t.hour, t.minute)
            data[d['name']].append(["", l.text, n.text, timestamp]) # first element is empty as we'll put the boat ID there
    return data

def write_csv(f, b, data):
    """
    Write the data to csv format
    :param f: file name to write csv to
    :param b: boat ids file. If it exists we'll read it otherwise we'll create it
    :param data: data to write
    :return:
    """

    boatids = {}
    maxid = 0
    if os.path.exists(b):
        with open(b, 'r') as bi:
            for r in csv.reader(bi):
                if len(r) > 0:
                    boatids[r[1]] = r[0]
                    if int(r[0]) > maxid:
                        maxid = int(r[0])

    with open(f, 'w') as fout:
        writer=csv.writer(fout)
        changed = False
        for boat in data:
            if boat not in boatids:
                maxid += 1
                boatids[boat] = maxid
                changed = True
            for d in data[boat]:
                d[0] = boatids[boat]
                writer.writerow(d)

    if changed:
        sys.stderr.write("Overwriting {} with new data!\n".format(b))
        with open(b, 'w') as bout:
            bwriter = csv.writer(bout)
            for b in boatids:
                bwriter.writerow([boatids[b], b])







if __name__ == '__main__':
    aparser = argparse.ArgumentParser(description="Convert an xml file to a csv file")
    aparser.add_argument('-f', help='XML file to convert', default="/home/redwards/Dropbox/Moving Around/PositionReports/sdpv16.xml")
    aparser.add_argument('-b', help='Boat ids file. if this doesnt exist we will make it', default="/home/redwards/Dropbox/Moving Around/PositionReports/sdpv16-boatids.txt")
    aparser.add_argument('-c', help='csv output file', default='/home/redwards/Desktop/temp.csv')
    aparser.add_argument('-v', help='verbose output', action="store_true")
    args = aparser.parse_args()

    data = read_xml(args.f)
    write_csv(args.c, args.b, data)
