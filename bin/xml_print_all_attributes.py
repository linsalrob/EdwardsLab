"""
Print all the attributes for all tags in an xml file
"""

import os
import sys
import argparse
import xml.etree.ElementTree as ET



def parse_child(child, tag_to_print, toprint):
    """
    A recursive method to just print tag and attribute labels
    :param child:
    :return:
    """

    if len(child.attrib) == 0:
        if toprint:
            print("{}".format(child.tag))

    if tag_to_print == child.tag or tag_to_print == 'all':
        toprint=True
    else:
        toprint = False

    # want to print attribute_name as this will become column headers

    for attr in child.attrib:
        if toprint:
            print("{}\t{}".format(child.tag, attr))
            #if 'Attribute' == child.tag and 'harmonized_name'==attr:
            #    print("{}\t{}\t{}".format(child.tag, attr, child.attrib['harmonized_name']))

    for gc in child:
        parse_child(gc, tag_to_print, toprint)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Print a redundant list of all the tags and their attributes (ie. use sort -u)')
    parser.add_argument('-f', help='filename for the xml file', required=True)
    parser.add_argument('-t', help='print a specific tag + its children (default = all tags)', default='all')
    args = parser.parse_args()


    toprint = True
    if args.t != 'all':
        toprint = False

    text=""


    with open(args.f, 'r') as fin:
        for l in fin:
            text = text + l
    # print(f"TEXT: {text}")
    biosampleset = ET.fromstring(text)
    for biosample in biosampleset:
        parse_child(biosample, args.t, toprint)
