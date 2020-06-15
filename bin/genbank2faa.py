"""
Convert a genbank file to a fasta amino acids file
"""

import os
import sys
import argparse
from Bio import SeqIO

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def genbank_to_faa(gbkf, complexheader=False, verbose=False):
    """
    Parse a genbank file
    :param gbkf:
    :param complexheader:
    :param verbose:
    :return:
    """

    for seq in SeqIO.parse(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                contine
            (start, stop, strand) = (feat.location.start.position, feat.location.end.position, feat.strand)
            prtmtd = {
                'EC_number': "",
                                'locus_tag': "",
                'note': "",
                'product': "",
                'protein_id': "",
                'ribosomal_slippage': "",
                'transl_table': 11,
                'translation': ""
            }

            cid = ">"
            if 'protein_id' in feat.qualifiers:
                cid += '|'.join(feat.qualifiers['protein_id'])
            elif 'locus_tag' in feat.qualifiers:
                    cid += "|".join(feat.qualifiers['locus_tag'])
            else:
                cid += seq.id + "." + str(feat.location)

            if complexheader:
                loc = f"{start}_{stop}"
                if strand < 0:
                    loc = f"{stop}_{start}"

                cid += f' [{seq.id}] [{seq.annotations["organism"]}] [{seq.id}_{loc}]'
                if 'product' in feat.qualifiers:
                    cid += f' [{feat.qualifiers["product"]}]'
                else:
                    cid += f' [hypothetical protein]'

            print(f"{cid}\n{feat.qualifiers["translation"]}")




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-g', help='genbank file', required=True)
    parser.add_argument('-c', help='complex identifier line', action='store_true')
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


