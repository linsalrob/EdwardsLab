"""

An assembly pipeline for Nanopore data using filtlong, miniasm, minipolish, and minimap.

Creates a fasta file of the assembly and generates some stats about it too.

To run on anthill use:
    mkdir sge_err sge_out
    snakemake  --snakefile nanopore_assembly.snakefile  --cluster 'qsub -q default -cwd -o sge_out -e sge_err -V'  -j 200 --latency-wait 60


Rob Edwards, Feb 2020

"""

import os
import re

READDIR = "fastq_pass"
OUTDIR = "assembly"

# a couple of REs used in filtering the files
barcodere = re.compile(r'^barcode\d+$')
fastqre   = re.compile(r'\.fastq$')

# define barcode to be only ^\w+$. This is important to stop cyclic 
# use of the wildcard as files are created.
wildcard_constraints:
   barcode = '\w+'


# filter only for those directories that have fastq files in them
def bfilter(wd):
    if not barcodere.search(wd):
        return False
    if next(filter(fastqre.search, os.listdir(os.path.join("fastq_pass", wd))), None):
        return True
    return False

# get a list of just the barcodes. Note that we currently ignore unclassified, though could add 
# if required by converting to a list and appending.
BARCODES = list(filter(bfilter, os.listdir(READDIR)))


def get_fastq_files(wildcards):
    """
    Create a list of the fastq files. It is important that we have a constraint
    on the wildcard here.
    Note that we only process fastq files in the output directory.
    """
    return list(map(lambda x: os.path.join(READDIR, wildcards.barcode, x), filter(fastqre.search, os.listdir(os.path.join(READDIR, wildcards.barcode)))))


rule all:
    """
    Note: For debugging purposes, I keep all the rule outputs, but 
    the last one is the currently active rule.

    I have experienced odd activity when specifying a specific rule
    vs. the rule all invocation
    """
    input:
        # expand(os.path.join(OUTDIR, "{barcode}.fastq.gz"), barcode=BARCODES)
        # expand(os.path.join(OUTDIR, "{barcode}.filtlong.fastq.gz"), barcode=BARCODES)
        # expand(os.path.join(OUTDIR, "{barcode}.overlaps.paf"), barcode=BARCODES)
        # expand(os.path.join(OUTDIR, "{barcode}.gfa"), barcode=BARCODES)
        # expand(os.path.join(OUTDIR, "{barcode}.polished.gfa"), barcode=BARCODES)
        # expand(expand(os.path.join(OUTDIR, "{barcode}.fasta"), barcode=BARCODES)
        expand(os.path.join(OUTDIR, "{barcode}.assembly.stats.tsv"), barcode=BARCODES)


rule concat:
    input:
        fqf = get_fastq_files
    output:
        o = os.path.join(OUTDIR, "{barcode}.fastq.gz")
    shell:
        "cat {input.fqf} | gzip >  {output.o}"


rule filtlong:
    input:
        os.path.join(OUTDIR, "{barcode}.fastq.gz")
    params:
        min_length = 2000,
        keep_percent = 90,
        target_bases = 3000000
    output:
        os.path.join(OUTDIR, "{barcode}.filtlong.fastq.gz")
    shell:
        "filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} --target_bases {params.target_bases} {input} | gzip > {output}"


rule minimap:
    input:
        os.path.join(OUTDIR, "{barcode}.filtlong.fastq.gz")
    output:
        os.path.join(OUTDIR, "{barcode}.overlaps.paf")
    threads: 1
    # threads: workflow.cores * 0.50
    shell:
        "minimap2 -t {threads} -x ava-ont {input} {input} > {output}"

rule miniasm:
    input:
        fq = os.path.join(OUTDIR, "{barcode}.filtlong.fastq.gz"),
        paf = os.path.join(OUTDIR, "{barcode}.overlaps.paf")
    output:
        os.path.join(OUTDIR, "{barcode}.gfa")
    shell:
        "miniasm -f {input.fq} {input.paf} > {output}"

rule minipolish:
    input:
        flfq = os.path.join(OUTDIR, "{barcode}.filtlong.fastq.gz"),
        gfa = os.path.join(OUTDIR, "{barcode}.gfa")
    output:
        os.path.join(OUTDIR, "{barcode}.polished.gfa")
    shell:
        'minipolish {input.flfq} {input.gfa} > {output}'

rule gfa2fasta:
    input:
        os.path.join(OUTDIR, "{barcode}.polished.gfa")
    output:
        os.path.join(OUTDIR, "{barcode}.fasta")
    params:
        b = '{barcode}'
    shell:
        "awk '/^S/{{print \">{params.b}_\"$2\"\\n\"$3}}' {input} > {output}"

rule assembly_stats:
    input:
        os.path.join(OUTDIR, "{barcode}.fasta")
    output:
        os.path.join(OUTDIR, "{barcode}.assembly.stats.tsv")
    shell:
        'python3 ~/bin/countfasta.py -t -f {input} > {output}'



