# A snakefile to process the results of a search SRA run. 
# Because I can't find the last one I wrote
# 
# Rob Edwards, May 2021

import os
import sys


TEMPDIR="/local/edwa0468"


rule all:
    input:
        "1",
        "empty_bams_removed",
        "sequence_hits.tsv"


rule extract_results:
    input:
        "results.zip"
    output:
        directory("1")
    shell:
        "unzip results.zip"

rule remove_empty:
    input:
        "1"
    output:
        temporary("empty_bams_removed")
    shell:
        """
        for F in $(find . -size -851c -name \*bam); do
            rm -f $F $F.bai;
        done;
        touch {output}
        """

rule idxstats:
    # we do this all in one rule so we still have access to the local
    # storage to remove the temporary files
    input:
        "1",
        "empty_bams_removed"
    output:
        "sequence_hits.tsv"
    conda:
        "envs/samtools.yaml"
    resources:
        mem_mb=640000
    shell:
        """
        ODIR=$(mktemp -d -p {TEMPDIR});
        C=0;
        for F in $(find . -name \*bam); do 
            SAMP=$(echo $F | sed -e 's/.\/[0-9]\+\///; s/.bam//');
            C=$((C+1));
            echo -e "Sample\t$SAMP" > $ODIR/out.$C;
            samtools idxstats $F | awk '$3 != 0' | grep -v ^\* | cut -f 1,3 | \
                    sed -e "s|^|$SAMP\t|" >> $ODIR/out.$C;
        done;
        joinlists.pl -z -h $ODIR/* > {output}
        rm -rf $ODIR
        """

