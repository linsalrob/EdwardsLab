"""
Convert a genbank file to a .ppt file

This is a simple genbank to tbl converter using BioPython
"""
import os
import sys
import binascii
import gzip
import re
from Bio import SeqIO

def is_gzip(gbkf):
    """
    Is the file compressed?
    :param gbkf:
    :return: true if compressed else false
    """

    with open(gbkf, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'


def feat_to_text(feat, qual):
    if qual in feat.qualifiers:
        return " ".join(feat.qualifiers[qual])
    return "-"

def convert_genbank(gbkf, printout=False, verbose=False):
    """
    Convert the genbank file to a table with the same columns of the ppt file
    :param gbkf: the genbank input file
    :param printout: print the table
    :param verbose: more output
    :return: the table
    """

    res = []

    if is_gzip(gbkf):
        handle = gzip.open(gbkf, 'rt')
    else:
        handle = open(gbkf, 'r')

    gire = re.compile('GI:(\d+)')
    cogre = re.compile('(COG\S+)')
    for seq in SeqIO.parse(handle, "genbank"):
        if verbose:
            sys.stderr.write(f"Parsing {seq.id}\n")
        for feat in seq.features:
            if feat.type != "CDS":
                continue

            gi = "-"
            if gire.match(feat_to_text(feat, 'db_xref')):
                gi = gire.match(feat_to_text(feat, 'db_xref'))[1]
            cog = "-"
            if cogre.match(feat_to_text(feat, 'product')):
                cog = cogre.match(feat_to_text(feat, 'product'))[1]

            if 'gene' in feat.qualifiers:
                gene = " ".join(feat.qualifiers)

            thisres = [
                f"{feat.location.start}..{feat.location.end}",
                "+" if feat.strand >= 0 else "-",
                (feat.length / 3) - 1,
                gi,
                feat_to_text(feat, 'gene'),
                feat_to_text(feat, 'locus_tag'),
                cog,
                feat_to_text(feat, 'product')
            ]

            if printout:
                print("\t".join(thisres))

            res.append(thisres)

    return res