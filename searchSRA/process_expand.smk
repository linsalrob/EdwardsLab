# A snakefile to process the results of a search SRA run. 
# Because I can't find the last one I wrote
# 
        
# Rob Edwards, May 2021

import os
import sys


# For more information on deciding which MAPQ value to set, please see
# http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
#MAPQ = 3
MAPQ = 10
RESULTS = "results"
FILTERED = f"filtered_Q{MAPQ}"
IDXSTATS = f"idxstats_Q{MAPQ}"
ABSTRACTS = "/home/edwa0468/GitHubs/EdwardsLab/searchSRA/searchSRA_abstracts.tsv.gz"

def expand_uncompressed_files(wildcards):
    checkpoint_output = checkpoints.extract_results.get(**wildcards).output[0]
    SMP, = glob_wildcards(os.path.join(checkpoint_output, "{sample}.bam"))
    file_names = expand(os.path.join(IDXSTATS, "{SAMPLE}.txt"),
                        SAMPLE=SMP, allow_missing=True)
    return file_names


def find_directories(wildcards):
    # checkpoint_output = checkpoints.extract_results.get(**wildcards).output[0]
    #filenames = expand_uncompressed_files(wildcards)
    idxstats_output = checkpoints.idx_stats_done.get(**wildcards).output.f
    DRS, SMPS, = glob_wildcards(os.path.join(idxstats_output, "{dir}", "{sample}.bam"))
    return expand(os.path.join(IDXSTATS, "{DIRS}"), DIRS=DRS)

rule all:
    input:
        "reads_per_project.tsv",
        "reads_per_sample.tsv"


checkpoint extract_results:
    input:
        "results.zip"
    output:
        directory(RESULTS)
    params:
        r = RESULTS
    shell:
        """
        mkdir -p {params.r};
        unzip results.zip -d {params.r};
        """

rule filter_reads:
    """ Filter reads by MAPQ score"""
    input:
        os.path.join(RESULTS, "{SAMPLE}.bam")
    output:
        os.path.join(FILTERED, "{SAMPLE}.bam")
    conda:
        "envs/samtools.yaml"
    params:
        results = RESULTS,
        filt = FILTERED,
        qual = MAPQ
    shell:
        """
        samtools view -bq {params.qual} -o {output} {input}
        samtools index {output}
        """


rule idxstats:
    input:
        os.path.join(FILTERED, "{SAMPLE}.bam")
    output:
        os.path.join(IDXSTATS, "{SAMPLE}.txt")
    conda:
        "envs/samtools.yaml"
    resources:
        mem_mb=6000
    params:
        filt = FILTERED,
        idx = IDXSTATS
    shell:
        """
        mkdir -p {params.idx}/$(dirname {wildcards.SAMPLE});
        samtools idxstats {input} | awk '$3 != 0' | \
            cut -f 1,3 | sed -e "s|^|{wildcards.SAMPLE}\t|" \
            > {output}
        """

checkpoint idx_stats_done:
    input:
        e = expand_uncompressed_files,
    output:
        d = IDXSTATS,
        f = temporary("idxstats.done")
    shell:
        "touch {output.f}"

rule merge:
    input:
        "idxstats.done",
        d = find_directories
    output:
        p = "reads_per_project.tsv",
        s = "reads_per_sample.tsv"
    params:
        abstracts = ABSTRACTS
    shell:
        """
        python3 \
                $HOME/GitHubs/EdwardsLab/searchSRA/merge_counts_abstracts.py \
                -d {input.d} -a {params.abstracts} -r {output.s} \
                -p {output.p}
        """
