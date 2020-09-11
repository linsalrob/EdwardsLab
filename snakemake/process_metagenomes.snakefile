"""

A metagenome assembly pipeline written by Rob Edwards, Jan 2020.

The idea is to start with pairs of reads and assemble those, and then 
combine all the assembled contigs and dereplicate using mmseqs, then 
map reads using bowtie2, and take the unmapped reads and reassemble.

Repeat that with the site data sets and then again with the sample
datasets.

A work in progress...

To run on the cluster use this command:

snakemake --configfile process_metagenomes.json -s ~/GitHubs/EdwardsLab/snakemake/process_metagenomes.snakefile --cluster 'qsub -cwd -o sge_out -e sge_err -V -pe make 8 -q default'  -j 200 --latency-wait 120

"""

import os
import sys
import socket


hostname = socket.gethostname()

sys.stderr.write(f"Running on: {hostname}\n")

configfile: 'process_metagenomes.json'


READDIR = config['directories']['Reads']
ASSDIR  = config['directories']['round1_assembly_output']
CRMDIR  = config['directories']['round1_contig_read_mapping']
UNASSM  = config['directories']['round2_unassembled_reads']
REASSM  = config['directories']['round2_assembly_output']
CCMO    = config['directories']['combined_contig_merging']
THREADS = config['threads']

if hostname.startswith('node'):
    hostname = 'anth'

if f"{hostname}_executables" not in config:
    sys.stderr.write(f"No executables defined for {hostname}_executables in process_metagenomes.json\n")
    sys.exit(-1)

config['executables'] = config[f"{hostname}_executables"]

# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(READDIR, '{sample}R1{extn}'))
if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)
if len(set(EXTENSIONS)) != 1:
    sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
    sys.stderr.write("\n\t".join(set(EXTENSIONS)))
    sys.stderr.write("\nWe don't know how to handle these\n")
    sys.exit(0)

FQEXTN = EXTENSIONS[0]
PATTERN_R1 = '{sample}R1' + FQEXTN
PATTERN_R2 = '{sample}R2' + FQEXTN

print(f"Samples are {SAMPLES}")

# this is the samples we have processed!

rule all:
    input:
        #os.path.join(ASSDIR, "round1_contigs.fa")
        #expand(os.path.join(CRMDIR, "{sample}.contigs.bam"), sample=SAMPLES)
        os.path.join(REASSM, "merged_contigs.fa"),
        os.path.join(CCMO, "flye.log"),


rule megahit_assemble:
    input:
        r1 = os.path.join(READDIR, PATTERN_R1),
        r2 = os.path.join(READDIR, PATTERN_R2)
    output:
        os.path.join(ASSDIR, "{sample}/final.contigs.fa"),
        os.path.join(ASSDIR, "{sample}/checkpoints.txt"),
        os.path.join(ASSDIR, "{sample}/done"),
        directory(os.path.join(ASSDIR, "{sample}/intermediate_contigs")),
        os.path.join(ASSDIR, "{sample}/log"),
        os.path.join(ASSDIR, "{sample}/options.json")
    params:
        odir = directory(os.path.join(ASSDIR, '{sample}'))
    threads: THREADS
    shell:
        'rmdir {params.odir}; {config[executables][assembler]} -1 {input.r1} -2 {input.r2} -o {params.odir} -t {threads}'

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
    threads: THREADS
    shell:
        # note that we build a large index by default, because
        # at some point we will end up doing that, and so this
        # always makes a large index
        'bowtie2-build --threads {threads} --large-index {input} {params.baseoutput}'
        

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
    threads: THREADS
    shell:
        'bowtie2 -x {params.contigs} -1 {params.r1} -2 {params.r2} --threads {threads} | samtools view -bh | samtools sort -o {output} -'

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
        os.path.join(REASSM, "checkpoints.txt"),
        os.path.join(REASSM, "done"),
        directory(os.path.join(REASSM, "intermediate_contigs")),
        os.path.join(REASSM, "log"),
        os.path.join(REASSM, "options.json")
    params:
        odir = os.path.join(REASSM)
    threads: THREADS
    shell:
        'rmdir {params.odir}; {config[executables][assembler]} -1 {input.r1} -2 {input.r2} -r {input.s0} -o {params.odir} -t {threads}'

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
        'python3 ~/bin/renumber_merge_fasta.py -f {input.new} -o {output.contigs} -i {output.ids} -v'

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
    threads: THREADS
    shell:
        """
        flye --meta --subassemblies {input.contigs} -o {CCMO} --threads {threads}
        """

