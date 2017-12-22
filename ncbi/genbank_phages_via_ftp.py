import StringIO
from ftplib import FTP
import gzip
from Bio import SeqIO

r = StringIO.StringIO()

def read_data(data):
    r.write(data)

ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
ftp.cwd('genbank/')
ftp.retrbinary('RETR gbphg3.seq.gz', r.write)

r.seek(0)

for seq in SeqIO.parse(gzip.GzipFile(fileobj=r), 'genbank'):
    print(seq.id + "\t" + seq.)