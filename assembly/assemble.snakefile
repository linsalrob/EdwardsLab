"""
Assemble the phage genomes.

This uses the minimap/miniasm/minipolish pipeline developed by Ryan Wick
https://github.com/rrwick/Minipolish

To run this on anthill use:
    snakemake --snakefile assemble.snakefile --cluster 'qsub -V -q default' --jobs 200 --latency-wait 60

[note that it is important to increase the latency wait to about 60 seconds. The default is 5 and the jobs will fail]

"""

from os.path import join
READDIR = "fastq_pass"

SAMPLES, = glob_wildcards(join(READDIR, '{sample}.fastq'))
# BARCODE, SAMPLES, = glob_wildcards(join(READDIR, '{barcode}', '{sample}.fastq'))
# print(f"BARCODE: {BARCODE}")
print(f"SAMPLES: {SAMPLES}")

rule all:
    input:
        expand(join(READDIR, '{sample}.assembly.stats.tsv'),  sample=SAMPLES)

rule minimap:
    input:
        '{sample}.fastq'
    output:
        '{sample}.overlaps.paf'
    threads: 1
    # threads: workflow.cores * 0.50
    shell:
        "minimap2 -t {threads} -x ava-ont {input} {input} > {output}"

rule miniasm:
    input:
        fq = '{sample}.fastq',
        paf = '{sample}.overlaps.paf'
    output:
        '{sample}.gfa'
    shell:
        "miniasm -f {input.fq} {input.paf} > {output}"

rule minipolish:
    input:
        fq = '{sample}.fastq',
        gfa = '{sample}.gfa'
    output:
        '{sample}.polished.gfa'
    threads: 1
    shell:
        'minipolish -t {threads} {input.fq} {input.gfa} > {output}'

rule gfa2fasta:
    input:
        '{sample}.polished.gfa'
    output:
        '{sample}.fasta'
    params:
        s = '{sample}'
    shell:
        "awk '/^S/{{print \">{params.s}_\"$2\"\\n\"$3}}' {input} > {output}"


rule assembly_stats:
    input:
        '{sample}.fasta'
    output:
        '{sample}.assembly.stats.tsv'
    shell:
        'python3 ~/bin/countfasta.py -t -f {input} > {output}'


