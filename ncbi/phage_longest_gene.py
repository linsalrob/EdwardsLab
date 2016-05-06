"""

Retrieve all the "longest" genes from phages in genbank ... what ever that means

"""

import os
import sys

import StringIO
from Bio import SeqIO
from ftplib import FTP


def process_a_file(filename):
    if 'phg' in filename:
        sys.stderr.write("Processing phage file: " + filename + "\n")
        print("Processing phage file: " + filename)
        r = StringIO
        ftp.retrbinary('RETR genbank/' + filename, r.write)
        for seq in SeqIO.read(r.getvalue(), 'genbank'):
            for feature in seq.features:
                if 'locus_tag' in feature.qualifiers:
                    lt = feature.qualifiers['locus_tag'][0]
                if 'translation' in feature.qualifiers:
                    print("{}\t{}".format(lt, len(feature.qualifiers['translation'][0])))
    else:
        sys.stderr.write("Skipped " + filename + "\n")

# first find all the phg files in genbank

ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
ftp.cwd('genbank/')
ftp.retrlines('NLST', process_a_file)
