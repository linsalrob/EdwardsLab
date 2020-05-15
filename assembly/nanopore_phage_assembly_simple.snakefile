"""

This is a simple version that only makes a fasta from Wick and a
polished flye fasta.

It does not make any stats

Rob Edwards, Feb 2020

"""

import os
import sys

## User defined options. 
# Yes, this should be in a config file, but its not

READDIR = "fastq"
OUTDIR  = "assembly"
FLYDIR  = "flye"


# this should be a path to your bacterial genomes that you want removed.
BACTERIALDNA = os.path.join(os.environ['HOME'], "phage/Sequencing/BacterialGenomes/Bacteria.fna")
if not os.path.exists(BACTERIALDNA):
    sys.stderr.write("FATAL: {BACTERIALDNA} does not exist\n")
    sys.exit()

# Some options for filtlong. You may want to change these
# https://github.com/rrwick/Filtlong
MIN_LENGTH = 2000,
KEEP_PERCENT = 90,
TARGET_BASES = 30000000

MINIPOLISH_ROUNDS = 10

# get a list of the fastq files
FASTQ, = glob_wildcards(os.path.join(READDIR, '{fastq}.fastq.gz'))


rule all:
    input:
        expand(os.path.join(OUTDIR, "{sample}_06.fasta"), sample=FASTQ),
        expand(os.path.join(OUTDIR, "{sample}_10_flye_polished.fasta"), sample=FASTQ),

"""
NOTE ABOUT THE OUTPUTS:

    The numbers are the step numbers. This enables natural sorting of
    the statistics at the end. Basically, everything is sorted by sample
    and then by number to be in the right order!

    Add a number on your output files.

    In addition, we add the number with _1 _2, etc so that we can
    use them later in the analysis without having to change their names.
"""


# first we want to remove any reads with "significant" hits to
# any of our genomes.

rule remove_host:
    input:
        os.path.join(READDIR, "{sample}.fastq.gz")
    output:
        os.path.join(OUTDIR, "{sample}_01.host_removed.fq.gz")
    shell:
        "minimap2 -ax map-ont {BACTERIALDNA} {input} | samtools fastq -f 4 | gzip > {output}"


rule filtlong:
    input:
        os.path.join(OUTDIR, "{sample}_01.host_removed.fq.gz")
    params:
        min_length = MIN_LENGTH, 
        keep_percent = KEEP_PERCENT,
        target_bases = TARGET_BASES
    output:
        os.path.join(OUTDIR, "{sample}_02.filtlong.fastq.gz")
    shell:
        "filtlong --min_length {params.min_length} --keep_percent {params.keep_percent} --target_bases {params.target_bases} {input} | gzip > {output}"


rule minimap:
    input:
        os.path.join(OUTDIR, "{sample}_02.filtlong.fastq.gz")
    output:
        os.path.join(OUTDIR, "{sample}_03.overlaps.paf")
    threads: 1
    # threads: workflow.cores * 0.50
    shell:
        "minimap2 -t {threads} -x ava-ont {input} {input} > {output}"

rule miniasm:
    input:
        fq = os.path.join(OUTDIR, "{sample}_02.filtlong.fastq.gz"),
        paf = os.path.join(OUTDIR, "{sample}_03.overlaps.paf")
    output:
        os.path.join(OUTDIR, "{sample}_04.miniasm.gfa")
    shell:
        "miniasm -f {input.fq} {input.paf} > {output}"

rule minipolish:
    input:
        # we run this with the original reads not the filtlong reads
        flfq = os.path.join(READDIR, "{sample}.fastq.gz"),
        #flfq = os.path.join(OUTDIR, "{sample}.2.filtlong.fastq.gz"),
        gfa = os.path.join(OUTDIR, "{sample}_04.miniasm.gfa")
    output:
        os.path.join(OUTDIR, "{sample}_05.polished.gfa")
    params:
        rds = MINIPOLISH_ROUNDS
    shell:
        'minipolish --rounds {params.rds} {input.flfq} {input.gfa} > {output}'

rule gfa2fasta:
    input:
        os.path.join(OUTDIR, "{sample}_05.polished.gfa")
    output:
        os.path.join(OUTDIR, "{sample}_06.fasta")
    params:
        b = '{sample}'
    shell:
        "awk '/^S/{{print \">{params.b}_\"$2\"\\n\"$3}}' {input} > {output}"



rule flye_assembly:
    input:
        os.path.join(OUTDIR, "{sample}_01.host_removed.fq.gz")
    output:
        os.path.join(FLYDIR, "{sample}", "assembly.fasta")
    params:
        odir = os.path.join(FLYDIR, "{sample}")
    shell:
        """
        flye --nano-raw {input} --asm-coverage 50 --genome-size 100k -o {params.odir} || touch {output}
        """


rule map_for_racon_polish_flye:
    """
    This was recommended by Ryan Wick. 

    First, we use minimap2 to map our original reads to the contigs, 
    and then we use racon to polish the contigs
    """
    input:
        reads = os.path.join(READDIR, "{sample}.fastq.gz"),
        conts = os.path.join(FLYDIR, "{sample}", "assembly.fasta")
    output:
        os.path.join(FLYDIR, "{sample}", "mapped_reads.paf")
    shell:
        """
        if [ -s {input.conts} ]; then
            minimap2 {input.conts} {input.reads} > {output}
        else
            touch {output}
        fi
        """

rule racon_polish_flye:
    input:
        reads = os.path.join(READDIR, "{sample}.fastq.gz"),
        conts = os.path.join(FLYDIR, "{sample}", "assembly.fasta"),
        maps  = os.path.join(FLYDIR, "{sample}", "mapped_reads.paf")
    output:
        os.path.join(FLYDIR, "{sample}", "polished_assembly.fasta")
    shell:
        """
        if [ -s {input.conts} ]; then
            racon {input.reads} {input.maps} {input.conts} > {output}
        else
            touch {output}
        fi
        """

rule move_flye_assemblies:
    input:
        p = os.path.join(FLYDIR, "{sample}", "polished_assembly.fasta")
    output:
        p = os.path.join(OUTDIR, "{sample}_10_flye_polished.fasta")
    shell:
        "cp {input.p} {output.p}"

