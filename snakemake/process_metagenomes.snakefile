"""

A metagenome assembly pipeline written by Rob Edwards, Jan 2020.

The idea is to start with pairs of reads and assemble those, and then 
combine all the assembled contigs and dereplicate using mmseqs, then 
map reads using bowtie2, and take the unmapped reads and reassemble.

Repeat that with the site data sets and then again with the sample
datasets.

To run on the cluster use this command:

snakemake --configfile process_metagenomes.json -s ~/GitHubs/EdwardsLab/snakemake/process_metagenomes.snakefile --profile sge

"""

import os
import sys
import socket




configfile: 'process_metagenomes.json'


READDIR = config['directories']['Reads']
ASSDIR  = config['directories']['round1_assembly_output']
CRMDIR  = config['directories']['round1_contig_read_mapping']
UNASSM  = config['directories']['round2_unassembled_reads']
REASSM  = config['directories']['round2_assembly_output']
CCMO    = config['directories']['combined_contig_merging']
PSEQDIR = config['directories']['prinseq']


# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(READDIR, '{sample}_R1.{extn}'))
if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)
if len(set(EXTENSIONS)) != 1:
    sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
    sys.stderr.write("\n\t".join(set(EXTENSIONS)))
    sys.stderr.write("\nWe don't know how to handle these\n")
    sys.exit(0)

FQEXTN = EXTENSIONS[0]
PATTERN_R1 = '{sample}_R1.' + FQEXTN
PATTERN_R2 = '{sample}_R2.' + FQEXTN


rule all:
    input:
        #os.path.join(ASSDIR, "round1_contigs.fa")
        #expand(os.path.join(CRMDIR, "{sample}.contigs.bam"), sample=SAMPLES)
        os.path.join(REASSM, "merged_contigs.fa"),
        os.path.join(CCMO, "flye.log"),


rule prinseq:
    input:
        r1 = os.path.join(READDIR, PATTERN_R1),
        r2 = os.path.join(READDIR, PATTERN_R2)
    output:
        r1 = os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq"),
        b1 = temporary(os.path.join(PSEQDIR, "{sample}_bad_out_R1.fastq")),
        b2 = temporary(os.path.join(PSEQDIR, "{sample}_bad_out_R2.fastq"))
    conda: "envs/prinseq.yaml"
    params:
        o = os.path.join(PSEQDIR, "{sample}")
    shell:
        """
            prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 \
                    -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
                    -ns_max_n 5  -trim_qual_type min -trim_qual_left 30 \
                    -trim_qual_right 30 -trim_qual_window 10 \
                    -threads {threads} \
                    -out_name {params.o} \
                    -fastq {input.r1} \
                    -fastq2 {input.r2};
        """

rule megahit_assemble:
    input:
        r1 = os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq")
    output:
        os.path.join(ASSDIR, "{sample}/final.contigs.fa"),
        temporary(os.path.join(ASSDIR, "{sample}/checkpoints.txt")),
        temporary(os.path.join(ASSDIR, "{sample}/done")),
        temporary(directory(os.path.join(ASSDIR, "{sample}/intermediate_contigs"))),
        temporary(os.path.join(ASSDIR, "{sample}/log")),
        temporary(os.path.join(ASSDIR, "{sample}/options.json"))
    params:
        odir = directory(os.path.join(ASSDIR, '{sample}'))
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/megahit.yaml"
    shell:
        """
        rmdir {params.odir}; 
        megahit -1 {input.r1} -2 {input.r2} -r {input.s1} -r {input.s2} \
                -o {params.odir} -t {resources.cpus}
        """

rule combine_contigs:
    """
    Here we take all contigs produced by megahit and combine them into a single (redundant)
    fasta file. This also creates a file called round1_contigs.ids that has the original
    contig ids
    """

    input:
        expand(os.path.join(ASSDIR, "{sample}/final.contigs.fa"), sample=SAMPLES)
    output:
        contigs = os.path.join(ASSDIR, "round1_contigs.fa"),
        ids = os.path.join(ASSDIR, "round1_contigs.ids")
    shell:
        'python3 ~/GitHubs/EdwardsLab/bin/renumber_merge_fasta.py -f {input} -o {output.contigs} -i {output.ids} -v'


rule index_contigs:
    """
    build the bowtie2 index of the contigs
    """
    input:
        os.path.join(ASSDIR, "round1_contigs.fa")
    params:
        baseoutput = os.path.join(CRMDIR, "round1_contigs")
    output:
        idx1 = temporary(os.path.join(CRMDIR, "round1_contigs.1.bt2l")),
        idx2 = temporary(os.path.join(CRMDIR, "round1_contigs.2.bt2l")),
        idx3 = temporary(os.path.join(CRMDIR, "round1_contigs.3.bt2l")),
        idx4 = temporary(os.path.join(CRMDIR, "round1_contigs.4.bt2l")),
        ridx1 = temporary(os.path.join(CRMDIR, "round1_contigs.rev.1.bt2l")),
        ridx2 = temporary(os.path.join(CRMDIR, "round1_contigs.rev.2.bt2l"))
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bowtie.yaml"
    shell:
        # note that we build a large index by default, because
        # at some point we will end up doing that, and so this
        # always makes a large index
        """
        mkdir -p {CRMDIR} && \
        bowtie2-build --threads {resources.cpus} --large-index {input} {params.baseoutput}
        """
        

rule map_reads:
    """
    Here we map the original reads back to the contigs
    """
    input:
        idx1 = os.path.join(CRMDIR, "round1_contigs.1.bt2l"),
        idx2 = os.path.join(CRMDIR, "round1_contigs.2.bt2l"),
        idx3 = os.path.join(CRMDIR, "round1_contigs.3.bt2l"),
        idx4 = os.path.join(CRMDIR, "round1_contigs.4.bt2l"),
        ridx1 = os.path.join(CRMDIR, "round1_contigs.rev.1.bt2l"),
        ridx2 = os.path.join(CRMDIR, "round1_contigs.rev.2.bt2l"),
        r1 = os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq")
    params:
        contigs = os.path.join(CRMDIR, "round1_contigs")
    output:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bowtie.yaml"
    shell:
        """
        bowtie2 -x {params.contigs} -1 {input.r1} -2 {input.r2} \
        -U {input.s1} -U {input.s2} --threads {resources.cpus} | \
        samtools view -bh | samtools sort -o {output} -
        """

"""
Using samtools to extract unmapped reads from the bam files
A samtools flag of -f 77 means the read is paird, neither the read nor mate is mapped, and the it is the first read in the pair, while a flag of -f 141 means the same except it is the second mate in the pair.
Then we use two flags: -f 4 (the read is unmapped) and -F 1 (the read is not paired) to find the single reads that are not mapped

See: https://edwards.sdsu.edu/research/command-line-deconseq/
"""

rule umapped_left_reads:
    """
    get left reads
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.R1.fastq")
    conda:
        "envs/bowtie.yaml"
    shell:
        """
        samtools view -h {input} | 
                awk 'BEGIN {{FS="\t"; OFS="\t"}} 
                {{if (/^@/ && substr($2, 3, 1)==":") {{print}} 
                else if (and($2, 0x1) && and($2, 0x40) && 
                (and($2, 0x4) || and($2, 0x8))) {{print}}}}' \
                | samtools bam2fq > {output}
        """

rule umapped_right_reads:
    """
    get right reads
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.R2.fastq")
    conda:
        "envs/bowtie.yaml"
    shell:
        """
        samtools view  -h {input} | 
                awk 'BEGIN {{FS="\t"; OFS="\t"}} 
                {{if (/^@/ && substr($2, 3, 1)==":") {{print}} 
                else if (and($2, 0x1) && and($2, 0x80) && 
                (and($2, 0x4) || and($2, 0x8))) {{print}}}}' \
                | samtools bam2fq > {output}
        """

rule umapped_single_reads:
    """
    get singletons
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.singles.fastq")
    conda:
        "envs/bowtie.yaml"
    shell:
            "samtools fastq -f 4 -F 1 {input} > {output}"

"""
We concatanate the unassembled reads into separate R1/R2/s files so
we can assmble them all together.

We also do this in 3 separate threads to take advantage of parallelization
and to make the command easier
"""

rule concatenate_R1_unassembled:
    """
    Start with R1 reads
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.R1.fastq"), sample=SAMPLES)
    output:
        os.path.join(UNASSM, "R1.unassembled.fastq")
    shell:
        "cat {input} > {output}"

rule concatenate_R2_unassembled:
    """
    Concat R2 reads
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.R2.fastq"), sample=SAMPLES)
    output:
        os.path.join(UNASSM, "R2.unassembled.fastq")
    shell:
        "cat {input} > {output}"

rule concatenate_single_unassembled:
    """
    Concate singletons
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.singles.fastq"),sample=SAMPLES)
    output:
        os.path.join(UNASSM, "single.unassembled.fastq")
    shell:
        "cat {input} > {output}"

rule assemble_unassembled:
    """
    assemble the unassembled reads
    """
    input:
        r1 = os.path.join(UNASSM, "R1.unassembled.fastq"),
        r2 = os.path.join(UNASSM, "R2.unassembled.fastq"),
        s0 = os.path.join(UNASSM, "single.unassembled.fastq")
    output:
        os.path.join(REASSM, "final.contigs.fa"),
        temporary(os.path.join(REASSM, "checkpoints.txt")),
        temporary(os.path.join(REASSM, "done")),
        temporary(directory(os.path.join(REASSM, "intermediate_contigs"))),
        temporary(os.path.join(REASSM, "log")),
        temporary(os.path.join(REASSM, "options.json"))
    params:
        odir = os.path.join(REASSM)
    conda:
        "envs/megahit.yaml"
    resources:
        mem_mb=64000,
        cpus=32
    shell:
        'rmdir {params.odir}; megahit -1 {input.r1} -2 {input.r2} -r {input.s0} -o {params.odir} -t {resources.cpus}'

"""
Combine all the contigs and use metaflye to merge the subassmblies
"""


rule concatentate_all_assemblies:
    """
    Again we take all contigs produced by megahit and combine them into a single (redundant)
    fasta file. This also creates a file called merged_contigs.ids that has the 
    contig ids. 
    """

    input:
        expand(os.path.join(ASSDIR, "{sample}/final.contigs.fa"), sample=SAMPLES),
        os.path.join(REASSM, "final.contigs.fa")
    output:
        contigs = os.path.join(REASSM, "merged_contigs.fa"),
        ids = os.path.join(REASSM, "merged_contigs.ids") 
    shell:
        'python3 ~/GitHubs/EdwardsLab/bin/renumber_merge_fasta.py -f {input} -o {output.contigs} -i {output.ids} -v'

rule merge_assemblies_with_flye:
    """
    Run flye on all the merged contigs fromt the last round to merge everything one more time
    """
    input:
        contigs = os.path.join(REASSM, "merged_contigs.fa")
    output:
        os.path.join(CCMO, "assembly.fasta"),
        os.path.join(CCMO, "assembly_graph.gfa"),
        os.path.join(CCMO, "assembly_graph.gv"),
        os.path.join(CCMO, "assembly_info.txt"),
        os.path.join(CCMO, "flye.log"),
    conda:
        "envs/flye.yaml"
    resources:
        mem_mb=64000,
        cpus=32
    shell:
        """
        flye --meta --subassemblies {input.contigs} -o {CCMO} --threads {resources.cpus}
        """

