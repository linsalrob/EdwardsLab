"""

Retrieve all the "longest" genes from phages in genbank ... what ever that means

"""

import os
import sys

import StringIO
from Bio import SeqIO
from ftplib import FTP
import gzip


r = StringIO.StringIO()

phgfiles = []

def read_data(data):
    r.write(data)

def phage_files(filename):
    if 'phg' in filename:
        phgfiles.append(filename)


def read_files():
    for filename in phgfiles:
        sys.stderr.write("Processing phage file: " + filename + "\n")
        ftp.retrbinary('RETR ' + filename, read_data)
        r.seek(0)

        for seq in SeqIO.parse(gzip.GzipFile(fileobj=r), 'genbank'):
            lt = ""
            for feature in seq.features:
                if 'protein_id' in feature.qualifiers:
                    lt = feature.qualifiers['protein_id'][0]
                if 'locus_tag' in feature.qualifiers:
                    lt = feature.qualifiers['locus_tag'][0]
                if 'translation' in feature.qualifiers:
                    print("{}\t{}\t{}".format(seq.id, lt, len(feature.qualifiers['translation'][0])))

    r.close()


# first find all the phg files in genbank

ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
ftp.cwd('genbank/')
ftp.retrlines('NLST', phage_files)
read_files()



