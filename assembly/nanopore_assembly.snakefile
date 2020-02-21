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
STATS = "stats"

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
        # expand(os.path.join(OUTDIR, "{barcode}.miniasm.gfa"), barcode=BARCODES)
        # expand(os.path.join(OUTDIR, "{barcode}.polished.gfa"), barcode=BARCODES)
        # expand(expand(os.path.join(OUTDIR, "{barcode}.fasta"), barcode=BARCODES)
        # expand(os.path.join(STATS, "{barcode}.assembly.stats.tsv"), barcode=BARCODES),
        os.path.join(STATS, "all_statistics.tsv")


"""
NOTE ABOUT THE OUTPUTS:

    The fasta/fastq files are named with the {barcode}.[\d+]

    The numbers are the step numbers. This enables natural sorting of
    the statistics at the end. Basically, everything is sorted by barcode
    and then by number to be in the right order!

    Add a number on your output files
"""


rule concat:
    input:
        fqf = get_fastq_files
    output:
        o = os.path.join(OUTDIR, "{barcode}.1.fastq.gz")
    shell:
        "cat {input.fqf} | gzip >  {output.o}"


rule filtlong:
    input:
        os.path.join(OUTDIR, "{barcode}.1.fastq.gz")
    params:
        min_length = 2000,
        keep_percent = 90,
        target_bases = 300000000
    output:
        os.path.join(OUTDIR, "{barcode}.2.filtlong.fastq.gz")
    shell:
        "filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} --target_bases {params.target_bases} {input} | gzip > {output}"


rule minimap:
    input:
        os.path.join(OUTDIR, "{barcode}.2.filtlong.fastq.gz")
    output:
        os.path.join(OUTDIR, "{barcode}.3.overlaps.paf")
    threads: 1
    # threads: workflow.cores * 0.50
    shell:
        "minimap2 -t {threads} -x ava-ont {input} {input} > {output}"

rule miniasm:
    input:
        fq = os.path.join(OUTDIR, "{barcode}.2.filtlong.fastq.gz"),
        paf = os.path.join(OUTDIR, "{barcode}.3.overlaps.paf")
    output:
        os.path.join(OUTDIR, "{barcode}.4.miniasm.gfa")
    shell:
        "miniasm -f {input.fq} {input.paf} > {output}"

rule minipolish:
    input:
        flfq = os.path.join(OUTDIR, "{barcode}.2.filtlong.fastq.gz"),
        gfa = os.path.join(OUTDIR, "{barcode}.4.miniasm.gfa")
    output:
        os.path.join(OUTDIR, "{barcode}.5.polished.gfa")
    shell:
        'minipolish {input.flfq} {input.gfa} > {output}'

rule gfa2fasta:
    input:
        os.path.join(OUTDIR, "{barcode}.5.polished.gfa")
    output:
        os.path.join(OUTDIR, "{barcode}.6.fasta")
    params:
        b = '{barcode}'
    shell:
        "awk '/^S/{{print \">{params.b}_\"$2\"\\n\"$3}}' {input} > {output}"

## Statistics on the sequences
rule count_concat_fastq:
    input:
        os.path.join(OUTDIR, "{barcode}.1.fastq.gz")
    output:
        os.path.join(STATS, "{barcode}_combined.tsv")
    shell:
        "python3 ~/EdwardsLab/bin/countfastq.py -t -f {input} > {output}"

rule count_filtlong:
    input:
        os.path.join(OUTDIR, "{barcode}.2.filtlong.fastq.gz")
    output:
        os.path.join(STATS, "{barcode}.filtlong.tsv")
    shell:
        "python3 ~/EdwardsLab/bin/countfastq.py -t -f {input} > {output}"


rule count_miniasm:
    input:
        os.path.join(OUTDIR, "{barcode}.4.miniasm.gfa")
    output:
        os.path.join(STATS, "{barcode}.miniasm.tsv")
    shell:
        "python3 ~/EdwardsLab/bin/countgfa.py -f {input} -t -v > {output}"

rule count_minipolish:
    input:
        os.path.join(OUTDIR, "{barcode}.5.polished.gfa")
    output:
        os.path.join(STATS, "{barcode}.polished.tsv")
    shell:
        "python3 ~/EdwardsLab/bin/countgfa.py -f {input} -t -v > {output}"

rule assembly_stats:
    input:
        os.path.join(OUTDIR, "{barcode}.6.fasta")
    output:
        os.path.join(STATS, "{barcode}.assembly.stats.tsv")
    shell:
        'python3 ~/bin/countfasta.py -t -f {input} > {output}'

rule combine_stats:
    input:
        expand(os.path.join(STATS, "{barcode}_combined.tsv"), barcode=BARCODES),
        expand(os.path.join(STATS, "{barcode}.filtlong.tsv"), barcode=BARCODES),
        expand(os.path.join(STATS, "{barcode}.miniasm.tsv"), barcode=BARCODES),
        expand(os.path.join(STATS, "{barcode}.polished.tsv"), barcode=BARCODES),
        expand(os.path.join(STATS, "{barcode}.assembly.stats.tsv"), barcode=BARCODES)
    output:
        os.path.join(STATS, "all_statistics.tsv")
    shell:
        "cat {input} >> {output}"

