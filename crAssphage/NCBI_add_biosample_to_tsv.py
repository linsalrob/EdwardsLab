"""
Read the biosample files and the biosample ID to the csv
"""

import os
import sys
import argparse
from roblib import bcolors

__author__ = 'Rob Edwards'


addnlheaders = ['alternate ID', 'bin_id', 'primer', 'sample', 'sampleid', 'sample_number', 'sample_processing', 'source', 'volunteer', 'yearmonth']

def read_biosids(bdfiles, verbose=False):
    """
    REad the biosample ID files.
    :param bdfiles: iterable of Files with biosample IDs
    :param verbose: more output
    :return: a dict of the sample_name and biosample ID
    """

    biosids = {}
    for fl in bdfiles:
        with open(fl, 'r') as f:
            for l in f:
                if l.startswith("Accession"):
                    continue
                p = l.rstrip('\n').split("\t")
                biosids[p[1]] = p[0]

    return biosids

def read_seqids(ifiles, biosids, verbose=False):
    """
    Read the ID files that links sample_name in the BioSamples
    to Sequence ID in the CSV files
    :param ifiles: the file of ids and sample_names. This is written by
    the -i flag of NCBI_submission.py
    :param biosids: biosample IDs
    :param verbose: more output
    :return: a dict of the seqid and the biosample ID
    """


    seqids = {}
    with open(ifiles, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if p[1] not in biosids:
                sys.stderr.write(f"{bcolors.RED}ERROR: |{p[1]}| does not have a biosample ID\n")
                continue
            seqids[p[0]] = biosids[p[1]]

    return seqids


def read_csv(csvindir, verbose=False):
    """
    Read the CSV files in csvindir and store some of the data

    :param csvindir: The directory with several csv input files
    :param verbose: more output
    """

    addnl = {}
    for inf in os.listdir(csvindir):
        with open(os.path.join(csvindir, inf), 'r') as f:
            headers = f.readline().rstrip('\n').split("\t")
            for l in f:
                p = l.rstrip('\n').split("\t")
                addnl[p[0]] = {}
                for t in addnlheaders:
                    i = headers.index(t)
                    addnl[p[0]][t] = p[i]
                addnl[p[0]]['sequence'] = p[headers.index('sequence')]
    return addnl

def add_to_everything(everythingf, addnl, seqids):
    """
    Read the everything file and add to it
    :param everythingf: the file from NCBI_submission.py
    :param addnl: the addtional data we missed
    :param seqids: the biosample IDs for the sequences
    """

    with open(everythingf, 'r') as f:
        for l in f:
            l = l.rstrip('\n')
            if l.startswith('sample_id'):
                ah = "\t".join(addnlheaders)
                print(f"BioProject_Accession\tBioSample_Accession\t{l}\t{ah}\tsequence")
                continue
            p = l.split("\t")
            if p[0] not in seqids:
                sys.stderr.write(f"{bcolors.RED}{p[0]} not in seqids{bcolors.ENDC}\n")
                continue
            if p[0] not in addnl:
                sys.stderr.write(f"{bcolors.RED}{p[0]} not in addnl{bcolors.ENDC}\n")
                continue

            ac = "\t".join([addnl[p[0]][h] for h in addnlheaders])
            print(f"PRJNA510571\t{seqids[p[0]]}\t{l}\t{ac}\t{addnl[p[0]]['sequence']}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-b', help='biosample files. You may specify more than one', required=True, action='append')
    parser.add_argument('-i', help='seqid -> seqname file from NCBI_submission.py', required=True)
    parser.add_argument('-c', help='csv input directory', required=True)
    parser.add_argument('-e', help='everything file from NCBI_submission.py', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    biosids = read_biosids(args.b, args.v)
    seqids  = read_seqids(args.i, biosids, args.v)
    addnl = read_csv(args.c, args.v)
    add_to_everything(args.e, addnl, seqids)

