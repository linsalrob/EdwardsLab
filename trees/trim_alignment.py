import argparse
import os
import sys

import re

__author__ = 'Rob Edwards'

"""
Trim an alignment file (in phylip format) based on the contents of the columns
"""


def parse_phylip_file(filename):
    """
    Parse the phylip alignment file

    :param filename:The file to parse
    :type filename: str
    :return: Two lists, one of the sequence ids and one of the alignments
    :rtype: list, list
    """

    seqids = []
    alignmentdata = []
    with open(filename, 'r') as f:
        l = f.readline()
        nalignments, alignlen = map(int, l.strip().split())
        currentline = 0
        for l in f:
            # sys.stderr.write("Current line: {} {}\n".format(currentline, l.strip()))
            # if currentline == nalignments - 1:
            #     currentline += 1
            #     continue
            if currentline == nalignments:
                currentline = 0
                continue
            if not l.startswith(' '):
                m = re.match('^\S+', l)
                if m:
                    seqids.append(m.group())
                else:
                    sys.stderr.write("Can't find a sequence ID at line {} : {}\n".format(currentline, l))
                l = re.sub('^\S+', '', l)
                alignmentdata.append("")
            alignmentdata[currentline] += l.strip().replace(' ', '')
            currentline += 1

    for a in alignmentdata:
        if len(a) != alignlen:
            sys.stderr.write("Alignment {} has length, but should be {}\n".format(len(a), alignlen))
    return seqids, alignmentdata



def print_alignment(seqids, alignment, linelength=65):
    """
    Print the alignment in phylip format

    :param seqids: The list of sequence ids
    :type seqids: list of str
    :param alignment: The list of alignment strings
    :type alignment: list of str
    :param linelength: the desired line length to print (default = 80)
    :type linelength: int
    :return:
    :rtype:
    """

    alignmentlen = len(alignment[0])
    naligns = len(alignment)

    lastend = 0
    print(" {} {}".format(naligns, alignmentlen))
    for i in range(len(seqids)):
        # print the sequence id and the correct number of spaces. max spaces = 11
        nspaces = 11-len(seqids[i])
        spaces = nspaces * ' '
        sys.stdout.write('{}{}'.format(seqids[i], spaces))
        end = 0
        for p in range(1, 6):
            start = end
            end = p * 10
            sys.stdout.write(alignment[i][start:end] + " ")
        sys.stdout.write("\n")
        lastend = end
    sys.stdout.write("\n")
    # now we just need to write out the rest of the data

    while lastend <= alignmentlen:
        for a in alignment:
            end = lastend
            sys.stdout.write(11 * ' ')
            for p in range(1, 6):
                start = end
                end = (p * 10) + lastend
                sys.stdout.write(a[start:end] + " ")
                # sys.stderr.write("{} to {}\n".format(start, end))
            sys.stdout.write("\n")
        sys.stdout.write("\n")
        lastend = end


def trim_alignment(alignment, cutoff):
    """
    Trim the alignment based on the number of informative residues

    :param alignment: The list of alignment strings
    :type alignment: list of str
    :param cutoff: The cutoff for informative residues
    :type cutoff: float
    :return: The revised list of alignments
    :rtype: list of str
    """

    alignmentlen = len(alignment[0])
    naligns = len(alignment)

    keepposn = []
    for i in range(alignmentlen):
        non_ir = 0 # non informative residues (i.e. '-')
        for a in alignment:
            if a[i] == '-':
                non_ir += 1
        if (1.0 * (naligns - non_ir) / naligns) >= cutoff:
            keepposn.append(i)


    newalignment = []
    for a in alignment:
        newalignment.append(''.join([a[i] for i in keepposn]))
    return newalignment





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trim an alignment file')
    parser.add_argument('-p', help='phylip alignment file', required=True)
    parser.add_argument('-c', help='minimum coverage as a fraction', required=True, type=float)
    args=parser.parse_args()

    seqids, alignment = parse_phylip_file(args.p)
    newalignment = trim_alignment(alignment, args.c)
    # newalignment = alignment
    print_alignment(seqids, newalignment)
