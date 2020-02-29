"""

A metagenome assembly pipeline written by Rob Edwards, Jan 2020.

The idea is to start with pairs of reads and assemble those, and then 
combine all the assembled contigs and dereplicate using mmseqs, then 
map reads using bowtie2, and take the unmapped reads and reassemble.

Repeat that with the site data sets and then again with the sample
datasets.

A work in progress...

"""

import os
import sys
import socket


hostname = socket.gethostname()

sys.stderr.write(f"Running on: {hostname}\n")

configfile: 'process_metagenomes.json'


wildcard_constraints:
   sample = '\w+'


READDIR = config['directories']['Reads']
ASSDIR  = config['directories']['round1_assembly_output']
TESTDIR = config['directories']['testdir']
CRMDIR  = config['directories']['round1_contig_read_mapping']
UNASSM  = config['directories']['round2_unassembled_reads']
REASSM  = config['directories']['round2_assembly_output']
RECRM   = config['directories']['round2_contig_read_mapping']

if hostname.startswith('node'):
    hostname = 'anth'

if f"{hostname}_executables" not in config:
    sys.stderr.write(f"No executables defined for {hostname}_executables in process_metagenomes.json\n")
    sys.exit(-1)

config['executables'] = config[f"{hostname}_executables"]

# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
SAMPLES, = glob_wildcards(os.path.join(READDIR, '{sample}_R1_001.fastq.gz'))
PATTERN_R1 = '{sample}_R1_001.fastq.gz'
PATTERN_R2 = '{sample}_R2_001.fastq.gz'

print(f"Samples are {SAMPLES}")

# this is the samples we have processed!

rule all:
    input:
        #os.path.join(ASSDIR, "round1_contigs.fa")
        #expand(os.path.join(CRMDIR, "{sample}.contigs.bam"), sample=SAMPLES)
        os.path.join(REASSM, "round2_contigs.fa"),
        expand(os.path.join(RECRM, "{sample}.contigs.bam"), sample=SAMPLES) 

rule megahit_assemble:
    input:
        r1 = os.path.join(READDIR, PATTERN_R1),
        r2 = os.path.join(READDIR, PATTERN_R2)
    output:
        os.path.join(ASSDIR, "{sample}/final.contigs.fa"),
        os.path.join(ASSDIR, "{sample}/checkpoints.txt"),
        os.path.join(ASSDIR, "{sample}/done"),
        os.path.join(ASSDIR, "{sample}/intermediate_contigs"),
        os.path.join(ASSDIR, "{sample}/log"),
        os.path.join(ASSDIR, "{sample}/options.json")
    params:
        odir = directory(os.path.join(ASSDIR, '{sample}')),
    shell:
        '{config[executables][assembler]} -1 {input.r1} -2 {input.r2} -o {params.odir}'

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
        'python3 ~/bin/renumber_merge_fasta.py -f {input} -o {output.contigs} -i {output.ids} -v'


rule index_contigs:
    """
    build the bowtie2 index of the contigs
    """
    input:
        os.path.join(ASSDIR, "round1_contigs.fa")
    params:
        baseoutput = os.path.join(CRMDIR, "round1_contigs")
    output:
        idx1 = os.path.join(CRMDIR, "round1_contigs.1.bt2l"),
        idx2 = os.path.join(CRMDIR, "round1_contigs.2.bt2l"),
        idx3 = os.path.join(CRMDIR, "round1_contigs.3.bt2l"),
        idx4 = os.path.join(CRMDIR, "round1_contigs.4.bt2l"),
        ridx1 = os.path.join(CRMDIR, "round1_contigs.rev.1.bt2l"),
        ridx2 = os.path.join(CRMDIR, "round1_contigs.rev.2.bt2l")
    shell:
        # note that we build a large index by default, because
        # at some point we will end up doing that, and so this
        # always makes a large index
        'bowtie2-build --large-index {input} {params.baseoutput}'
        

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
        ridx2 = os.path.join(CRMDIR, "round1_contigs.rev.2.bt2l")
    params:
        r1 = os.path.join(READDIR, PATTERN_R1),
        r2 = os.path.join(READDIR, PATTERN_R2),
        contigs = os.path.join(CRMDIR, "round1_contigs")
    output:
        os.path.join(CRMDIR, "{sample}.contigs.bam") 
    shell:
        'bowtie2 -x {params.contigs} -1 {params.r1} -2 {params.r2} | samtools view -bh | samtools sort -o {output} -'

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
    shell:
            "samtools fastq -f 4 -F 1 {input} > {output}"


rule assemble_unassembled:
    """
    assemble the unassembled reads
    """
    input:
        r1 = os.path.join(UNASSM, "{sample}.unassembled.R1.fastq"),
        r2 = os.path.join(UNASSM, "{sample}.unassembled.R2.fastq"),
        s0 = os.path.join(UNASSM, "{sample}.unassembled.singles.fastq")
    output:
        os.path.join(REASSM, "{sample}/final.contigs.fa"),
        os.path.join(REASSM, "{sample}/checkpoints.txt"),
        os.path.join(REASSM, "{sample}/done"),
        os.path.join(REASSM, "{sample}/intermediate_contigs"),
        os.path.join(REASSM, "{sample}/log"),
        os.path.join(REASSM, "{sample}/options.json")
    params:
        odir = os.path.join(REASSM, "{sample}")
    shell:
        '{config[executables][assembler]} -1 {input.r1} -2 {input.r2} -r {input.s0} -o {params.odir}'

"""
Start round 2.

Now we are going to repeat combining the contigs and mapping the sequences again
"""


rule concatentate_all_assemblies:
    """
    Again we take all contigs produced by megahit and combine them into a single (redundant)
    fasta file. This also creates a file called round2_contigs.ids that has the 
    contig ids. 
    """

    input:
        new = expand(os.path.join(REASSM, "{sample}/final.contigs.fa"), sample=SAMPLES)
    output:
        contigs = os.path.join(REASSM, "round2_contigs.fa"),
        ids = os.path.join(REASSM, "round2_contigs.ids") 
    shell:
        'python3 ~/bin/renumber_merge_fasta.py -f {input.ori} {input.new} -o {output.contigs} -i {output.ids} -v'


rule index_contigs_round2:
    """
    build the bowtie2 index of the contigs
    """
    input:
        os.path.join(REASSM, "round2_contigs.fa")
    params:
        baseoutput = os.path.join(RECRM, "round2_contigs")
    output:
        idx1 = os.path.join(RECRM, "round2_contigs.1.bt2l"),
        idx2 = os.path.join(RECRM, "round2_contigs.2.bt2l"),
        idx3 = os.path.join(RECRM, "round2_contigs.3.bt2l"),
        idx4 = os.path.join(RECRM, "round2_contigs.4.bt2l"),
        ridx1 = os.path.join(RECRM, "round2_contigs.rev.1.bt2l"),
        ridx2 = os.path.join(RECRM, "round2_contigs.rev.2.bt2l")
    shell:
        # note that we build a large index by default, because
        # at some point we will end up doing that, and so this
        # always makes a large index
        'bowtie2-build --large-index {input} {params.baseoutput}'
        

rule map_reads_round2:
    """
    Here we map the original reads back to the contigs
    """
    input:
        idx1 = os.path.join(RECRM, "round2_contigs.1.bt2l"),
        idx2 = os.path.join(RECRM, "round2_contigs.2.bt2l"),
        idx3 = os.path.join(RECRM, "round2_contigs.3.bt2l"),
        idx4 = os.path.join(RECRM, "round2_contigs.4.bt2l"),
        ridx1 = os.path.join(RECRM, "round2_contigs.rev.1.bt2l"),
        ridx2 = os.path.join(RECRM, "round2_contigs.rev.2.bt2l")
    params:
        r1 = os.path.join(READDIR, PATTERN_R1),
        r2 = os.path.join(READDIR, PATTERN_R2),
        contigs = os.path.join(RECRM, "round2_contigs")
    output:
        os.path.join(RECRM, "{sample}.contigs.bam") 
    shell:
        'bowtie2 -x {params.contigs} -1 {params.r1} -2 {params.r2} | samtools view -bh | samtools sort -o {output} -'



rule test:
    """
    Just test the snakefile
    """

    shell:
        "echo hello world!"

