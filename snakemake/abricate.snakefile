"""
Run [abricate](https://github.com/tseemann/abricate) with all 
options on a directory of sequence files
"""

import os
import sys
import subprocess

configfile: "abricate.yaml"


indir = config['seq_files']
outdir = config['abricate_output']

# read the known databases from abricate

proc = subprocess.Popen(['abricate','--list'],stdout=subprocess.PIPE, encoding='utf-8')
databases = list(filter(lambda x: x != "DATABASE", [p[0] for p in [l.strip().split("\t") for l in proc.stdout]]))

SEQS, = glob_wildcards(os.path.join(indir, '{seq}'))


rule all:
    input:
        expand(os.path.join(outdir, "{sample}.{db}.abricate.tsv"), sample=SEQS, db=databases)

rule abricate:
    input:
        os.path.join(indir, "{sample}")
    output:
        os.path.join(outdir, "{sample}.{db}.abricate.tsv")
    params:
        db = "{db}"
    shell:
        "abricate --noheader --nopath  --db {params.db} {input} > {output}"
