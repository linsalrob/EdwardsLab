import os
import sys

files = {}
for f in os.listdir('fastq'):
    sid = f.split('_')[0]
    if sid not in files:
        files[sid]=set()
    files[sid].add(f)

for s in files:
    if len(files[s]) == 1:
        print("fastq/" + files[s].pop())
    else:
        print("fastq/" + " fastq/".join(files[s]))
