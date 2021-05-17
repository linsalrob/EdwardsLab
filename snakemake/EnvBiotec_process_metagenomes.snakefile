"""

This is a snakefile for the processes that are done for the 2021 Environmental Biotech Class for Jess


Rob Edwards, January 2021

"""


import os
import sys


# set this to whatever the name of your directory
# with the reads is. If you are following along with the
# tutorial, you can leave this as fastq
READDIR = 'processed_by_rob/'

SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(READDIR, '{sample}.{extensions}'))

# just get the first file extension as we don't need to iterate all of them
file_extension = EXTENSIONS[0]


# just check there is something to actually do!
if len(SAMPLES) == 0:
    sys.stderr.write("FATAL: We could not detect any samples at all.\n")
    sys.stderr.write(f"Do you have a directory called {READDIR} with some fastq files in it?\n")
    sys.exit()


rule all:
    input:
        expand(
            [
                os.path.join("results", "{sample}", "prinseq", "{sample}.good.fastq"), 
                os.path.join("results", "{sample}", "focus", "output_All_levels.csv"),
                os.path.join("results", "{sample}", "superfocus", "output_all_levels_and_function.xls"),
                os.path.join("results", "{sample}", "kraken", "{sample}.report.tsv"),
                os.path.join("results", "{sample}", "kraken", "{sample}.output.tsv"),
            ],
               sample=SAMPLES),
        os.path.join("results", "all_taxonomy.tsv")



rule run_prinseq:
    input:
        r1 = os.path.join(READDIR, "{sample}." + file_extension),
    output:
        d = directory(os.path.join("results", "{sample}", "prinseq")),
        r1 = os.path.join("results", "{sample}", "prinseq", "{sample}.good.fastq"),
        b1 = temporary(os.path.join("results", "{sample}", "prinseq", "{sample}.bad.fastq")),
    resources:
        cpus=2,
        mem_mb=4000
    shell:
        """
        mkdir -p {output.d} &&
        prinseq++ -min_len 60 -min_qual_mean 1 -ns_max_n 1 -derep 1 \
            -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
            -out_good {output.r1} -out_bad {output.b1} \
            -fastq {input.r1} \
            -threads {resources.cpus}
        """

rule cp_psq_good:
    input:
        r1 = os.path.join("results", "{sample}", "prinseq", "{sample}.good.fastq")
    output:
        d = directory(os.path.join("results", "{sample}", "prinseq_good")),
        f = os.path.join("results", "{sample}", "prinseq_good", "{sample}.good.fastq")
    params:
        d = os.path.join("results", "{sample}", "prinseq_good")
    shell:
        """
        mkdir -p {params.d} && cp {input.r1} {output.f}
        """

rule run_focus:
    input:
        os.path.join("results", "{sample}", "prinseq_good")
    output:
        d = directory(os.path.join("results", "{sample}", "focus")),
        a = os.path.join("results", "{sample}", "focus", "output_All_levels.csv"),
        f = temporary(os.path.join("results", "{sample}", "focus", "output_Family_tabular.csv")),
        k = temporary(os.path.join("results", "{sample}", "focus", "output_Kingdom_tabular.csv")),
        p = temporary(os.path.join("results", "{sample}", "focus", "output_Phylum_tabular.csv")),
        s = temporary(os.path.join("results", "{sample}", "focus", "output_Strain_tabular.csv")),
        c = temporary(os.path.join("results", "{sample}", "focus", "output_Class_tabular.csv")),
        g = temporary(os.path.join("results", "{sample}", "focus", "output_Genus_tabular.csv")),
        o = temporary(os.path.join("results", "{sample}", "focus", "output_Order_tabular.csv")),
        sp = temporary(os.path.join("results", "{sample}", "focus", "output_Species_tabular.csv"))
    resources:
        cpus=2,
        mem_mb=4000
    shell:
        """
        focus -q {input} -o {output.d} -t {resources.cpus}
        """



rule run_superfocus:
    input:
        os.path.join("results", "{sample}", "prinseq_good")
    output:
        d = directory(os.path.join("results", "{sample}", "superfocus")),
        a = os.path.join("results", "{sample}", "superfocus", "output_all_levels_and_function.xls"),
        m8 = os.path.join("results", "{sample}", "superfocus", "{sample}.good.fastq_alignments.m8"),
        l1 = temporary(os.path.join("results", "{sample}", "superfocus", "output_subsystem_level_1.xls")),
        l2 = temporary(os.path.join("results", "{sample}", "superfocus", "output_subsystem_level_2.xls")),
        l3 = temporary(os.path.join("results", "{sample}", "superfocus", "output_subsystem_level_3.xls")),
    resources:
        cpus=4,
        mem_mb=16000
    shell:
        """
        superfocus -q {input} -dir {output.d} -a diamond -t {resources.cpus}
        """

rule run_kraken:
    input:
        r1 = os.path.join("results", "{sample}", "prinseq", "{sample}.good.fastq"),
    output:
        rt = os.path.join("results", "{sample}", "kraken", "{sample}.report.tsv"),
        ot = os.path.join("results", "{sample}", "kraken", "{sample}.output.tsv")
    resources:
        cpus=4,
        mem_mb=400000
    shell:
        """
        kraken2 --report {output.rt} \
                --output {output.ot} \
                --threads {resources.cpus} \
                {input.r1}
        """

rule superfocus_taxonomy:
    input:
        m8 = os.path.join("results", "{sample}", "superfocus", "{sample}.good.fastq_alignments.m8"),
    output:
        os.path.join("results", "{sample}", "superfocus", "{sample}_good_out.taxonomy")
    params:
        db = "/home/edwa0468/ncbi/taxonomy/taxonomy.sqlite3"
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/superfocus_all/superfocus_to_taxonomy.py -f {input} --tophit -d {params.db} > {output}
        """

rule count_sf_taxonomy:
    input:
        os.path.join("results", "{sample}", "superfocus", "{sample}_good_out.taxonomy")
    output:
        os.path.join("results", "{sample}", "superfocus", "{sample}_tax_counts.tsv")
    params:
        s = "{sample}"
    shell:
        """
        cut -d$'\t' -f 2- {input} | sed -e 's/\t/|/g' | \
        awk -v F={params.s} 'BEGIN {{print "Superkingdom|Phylum|Class|Order|Family|Genus|Species|Taxid\t"F}} s[$0]++ {{}} END {{ for (i in s) print i"\t"s[i] }}' \
        > {output}
        """

rule join_superfocus_taxonomy:
    input:
        expand(os.path.join("results", "{smps}", "superfocus", "{smps}_tax_counts.tsv"), smps=SAMPLES)
    output:
        os.path.join("results", "all_taxonomy.tsv")
    shell:
        """
        ~/GitHubs/EdwardsLab/bin/joinlists.pl -h -z {input} > {output}
        """



