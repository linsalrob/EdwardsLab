"""
Run virsorter on a genbank file

"""


import os
import sys


GBDIR = "genbank"
FADIR = "fasta"

SAMPLES, =  glob_wildcards(os.path.join(GBDIR, '{sample}.gbf'))



rule all:
    input:
        expand(os.path.join(FADIR, "{sample}.fna"), sample=SAMPLES)


rule genbank2fasta:
    input:
        os.path.join(GBDIR, "{sample}.gbf")
    output:
        os.path.join(FADIR, "{sample}.fna")
    shell:
        "any2fasta {input} > {output}"


rule run_virsorter:
    input:
        os.path.join(FADIR, "{sample}.fna")
    conda:
        "virsorter"
    output:
        os.path.join(VIRDIR, "{sample}", 




