import os
import subprocess
import sys

__author__ = 'Rob Edwards'

print('CWD'  + os.getcwd())

fname =  os.path.join(os.environ['HOME'], 'Dropbox/Metagenomics/51.hits_small.fa.lrz')

f = subprocess.Popen(['/usr/bin/lrunzip', '-q', '-d', '-f', '-o-', fname], stdout=subprocess.PIPE).stdout

for l in f:
    print("READ: {}".format(l))