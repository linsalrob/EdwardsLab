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

def read_data(data):
    r.write(data)

def process_a_file(filename):
    if 'phg' in filename:
        sys.stderr.write("Processing phage file: " + filename + "\n")
        ftp.retrbinary('RETR genbank/' + filename, read_data)
        r.seek(0)

        for seq in SeqIO.parse(gzip.GzipFile(fileobj=r), 'genbank'):
            for feature in seq.features:
                if 'locus_tag' in feature.qualifiers:
                    lt = feature.qualifiers['locus_tag'][0]
                if 'translation' in feature.qualifiers:
                    print("{}\t{}".format(lt, len(feature.qualifiers['translation'][0])))
    else:
        sys.stderr.write("Skipped " + filename + "\n")

    r.close()

# first find all the phg files in genbank

ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
ftp.cwd('genbank/')
ftp.retrlines('NLST', process_a_file)
