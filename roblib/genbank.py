"""
Read a genbank file and do things with it!
"""

import binascii
import gzip
from Bio import SeqIO

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def is_gzip(gbkf):
    """
    Is the file compressed?
    :param gbkf:
    :return: true if compressed else false
    """

    with open(gbkf, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'


def genbank_seqio(gbkf, verbose=False):
    """
    Get the parser stream
    :param gbkf: genbank file
    :param verbose:
    :return:
    """

    if is_gzip(gbkf):
        handle = gzip.open(gbkf, 'rt')
    else:
        handle = open(gbkf, 'r')
    return SeqIO.parse(handle, "genbank")


def genbank_to_faa(gbkf, complexheader=False, verbose=False):
    """
    Parse a genbank file
    :param gbkf:
    :param complexheader:
    :param verbose:
    :return: a dict of the sequences
    """

    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
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

                cid += f' [{seq.id}] '
                if 'organism' in seq.annotations:
                    cid += f' [{seq.annotations["organism"]}]'
                cid += f' [{seq.id}_{loc}]'
                if 'product' in feat.qualifiers:
                    cid += f' {feat.qualifiers["product"][0]}'
                else:
                    cid += f' [hypothetical protein]'

            yield cid, feat.qualifiers['translation'][0]


def genbank_to_orfs(gbkf, complexheader=False, verbose=False):
    """
    Parse a genbank file
    :param gbkf:
    :param complexheader:
    :param verbose:
    :return: a dict of the sequences
    """

    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
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

                cid += f' [{seq.id}] '
                if 'organism' in seq.annotations:
                    cid += f' [{seq.annotations["organism"]}]'
                cid += f' [{seq.id}_{loc}]'
                if 'product' in feat.qualifiers:
                    cid += f' {feat.qualifiers["product"][0]}'
                else:
                    cid += f' [hypothetical protein]'

            yield cid, feat.seq




def genbank_to_fna(gbkf):
    """
    Parse a genbank file
    :param gbkf: genbank file
    :param verbose:
    :return: a dict of the sequences
    """

    for seq in genbank_seqio(gbkf):
        yield seq.id, seq.seq
