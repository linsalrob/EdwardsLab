"""

Retrieve all the "longest" genes from phages in genbank ... what ever that means

"""

import os
import sys

import StringIO
from ftplib import FTP

# first find all the phg files in genbank

ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
ftp.cwd('genbank/')
for f in ftp.retrlines('LIST'):
    if 'phg' in f:
        print("Phage file: " + f)