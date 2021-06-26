# A snakefile to process the results of a search SRA run. 
# Because I can't find the last one I wrote
# 
# Rob Edwards, May 2021

import os
import sys


TEMPDIR = "/local/edwa0468"
RESULTS = "results"
FILTERED = "filtered"
MAPQ = 2
IDXSTATS = "idxstats"
ABSTRACTS = "/home/edwa0468/GitHubs/EdwardsLab/searchSRA/searchSRA_abstracts.tsv.gz"

rule all:
    input:
        "reads_per_project.tsv",
        "reads_per_sample.tsv"


rule extract_results:
    input:
        "results.zip"
    output:
        directory(os.path.join(RESULTS, "1"))
    params:
        results = RESULTS
    shell:
        """
        mkdir -p {params.results} && 
        unzip results.zip -d {params.results}
        """

rule remove_empty:
    input:
        os.path.join(RESULTS, "1")
    output:
        temporary("empty_bams_removed")
    params:
        results = RESULTS
    shell:
        """
        for F in $(find {params.results} -size -900c -name \*bam); do
            rm -f $F $F.bai;
        done;
        touch {output}
        """

rule filter_reads:
    """ Filter reads by MAPQ score"""
    input:
        os.path.join(RESULTS, "1")
    output:
        directory(FILTERED)
    conda:
        "envs/samtools.yaml"
    params:
        results = RESULTS,
        filt = FILTERED,
        qual = MAPQ
    shell:
        """
        mkdir -p {params.filt};
        for BAMFILE in $(find {params.results} -name \*bam); do
            BAM=$(basename $BAMFILE);
            samtools view -bq {params.qual} -o {params.filt}/$BAM $BAMFILE;
            samtools index {params.filt}/$BAM;
        done
        """


rule idxstats:
    input:
        FILTERED
    output:
        IDXSTATS
    conda:
        "envs/samtools.yaml"
    resources:
        mem_mb=640000
    params:
        filt = FILTERED,
        idx = IDXSTATS
    shell:
        """
        C=0;
        mkdir -p {params.idx};
        for BAMFILE in $(find {params.filt} -name \*bam); do
            FNAME=$(basename $BAMFILE);
            SAMP=$(echo $FNAME | sed -e 's/.bam//');
            C=$((C+1));
            samtools idxstats $BAMFILE | awk '$3 != 0' | \
                    grep -v ^\* | cut -f 1,3 | sed -e "s|^|$SAMP\t|" \
                    >> {params.idx}/out.$C;
        done;
    """
    
rule merge:
    input:
        IDXSTATS
    output:
        p = "reads_per_project.tsv",
        s = "reads_per_sample.tsv"
    params:
        abstracts = ABSTRACTS
    shell:
        """
        "python3 \
                $HOME/GitHubs/EdwardsLab/searchSRA/merge_counts_abstracts.py \
                -d {input} -a {params.abstracts} -r {output.s} \
                -p {output.p}
        """
