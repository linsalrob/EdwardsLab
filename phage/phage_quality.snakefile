"""
A pipeline to check the quality of phage contigs.

We expect:
    1. An input directory of contigs
    2. Within that directory each fasta file is
       a single phage genome. You may have as
       many contigs per fasta file as you wish

We will:
    1. Use phanotate to call ORFs
    2. Compare those ORFs to all phage proteins
        i.    compare query length and subject length
        ii.   look for adjacent genes that map to the same
              protein
    3. Compare those ORFs to bacterial proteins (probably a small subset ... TBD)
    4. Run viromeQC on the contigs
    5. Measure a bunch of stats about your contigs
        i.    N50 (add L50/U50?)
        ii.   N75
        iii.  How many frameshifted proteins
        iv.   genome length too large
        v.    genome length too small
        vi.   abnormal gene to sequence ratio
        vii.  low gene count
        viii. average gene length
        ix.   compare using translated blast (blastp or diamond??) and check for same hits on
              different strands. This should be similar to (iii)
"""

"""
Requirements:
    1. phage_quality_config.yaml -- the configuration file
    2. phanotate.py
    3. transeq from EMBOSS
    4. blastp
    5. blastx (not yet required ... may use diamond)

To run this on anthill start with:

    snakemake --configfile config.yaml -s ~/GitHubs/EdwardsLab/phage/phage_quality.snakefile --cluster 'qsub -cwd -o sge_out -e sge_err -V -pe make 8 ' --latency-wait 60 -j 600



"""


"""
TODO
- we need to provide a phages.faa database and format it for blast
- we need a bacterial database (nr? amphora?)
"""

import os
import sys
from roblib import bcolors, stream_fasta, median

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()

CONTIGS = config['paths']['contigs']
ORFS    = config['paths']['orfs']
BLAST   = config['paths']['blast']
STATS   = config['paths']['statistics']

PHAGEDB = os.path.join(config['paths']['databases'], config['databases']['phage_proteins'])
if not os.path.exists(PHAGEDB):
    sys.stderr.write(f"FATAL: {PHAGEDB} not found. Please check your paths\n")
    sys.exit()


SAMPLES, = glob_wildcards(os.path.join(CONTIGS, '{sample}.fasta'))


def count_adjacent_orfs(sample, fastafile, blastfile, adjacentout, nohitsout, searchtype):
    """
    Count hits where two adjacent orfs match
    to the same protein
    """
    sys.stderr.write(f"{bcolors.GREEN}Counting adjacent ORFs for {sample} and {searchtype}{bcolors.ENDC}\n")
    orfs = []
    with open(fastafile, 'r') as f:
        for l in f:
            if l.startswith('>'):
                p = l.split()
                seqid = p[0].replace('>', '', 1)
                orfs.append(seqid)
    hits = {}
    with open(blastfile, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if p[0] not in hits:
                hits[p[0]] = set()
            hits[p[0]].add(p[1])
    
    adjacent = 0
    nohits   = 0
    for i,o in enumerate(orfs):
        if i >= len(orfs)-1:
            continue
        if o not in hits:
            nohits += 1 
            continue
        if orfs[i+1] not in hits:
            continue
        for h in hits[o]:
            if h in hits[orfs[i+1]]:
                adjacent += 1
                break
    with open(adjacentout, 'w') as out:
        out.write(f"{sample}\tNumber of orfs with adjacent {searchtype} similarities\t")
        out.write(f"[adjacent orfs, total orfs, fraction adjacent]\t")
        out.write(f"{adjacent}\t{len(orfs)}\t{adjacent/len(orfs)}\n")
    with open(nohitsout, 'w') as out:
        out.write(f"{sample}\tNumber of orfs with no hits\t")
        out.write(f"[nohits, total orfs, fraction no hits]\t")
        out.write(f"{nohits}\t{len(orfs)}\t{nohits/len(orfs)}\n")


def av_protein_lengths(sample, blastfile, fractionout, summaryout, searchtype):
    """
    Calculate the average length of the best hit of all the proteins
    """
    q = {}
    av = []
    sys.stderr.write(f"{bcolors.GREEN}Average protein lengths for {sample} and {searchtype}{bcolors.ENDC}\n")

    with open(blastfile, 'r') as f:
        with open(fractionout, 'w') as out:
            for l in f:
                p = l.strip().split("\t")
                if p[0] in q:
                    continue
                q[p[0]] = int(p[12])/int(p[13])
                av.append(q[p[0]])
                out.write(f"{p[0]}\t{q[p[0]]}\n")
    with open(summaryout, 'w') as out:
        out.write(f"{sample}\tAverage {searchtype} protein lengths\t")
        out.write("[num orfs, median proportional length, averagh proportional length]\t")
        out.write(f"{len(av)}\t{median(av)}\t{sum(av)/len(av)}\n")

        
def coding_versus_noncoding(sample, contigs, orfs, outputfile):
    """
    Count the number of coding vs. noncoding bases
    """
    sys.stderr.write(f"{bcolors.GREEN}Counting coding vs. non coding bases{bcolors.ENDC}\n")

    # read the DNA sequences
    dnalen = 0
    for seqid, seq in stream_fasta(contigs):
        dnalen += len(seq)

    # read the protein sequences
    protlen = 0
    for seqid, seq in stream_fasta(orfs):
        protlen += len(seq)
    
    with open(outputfile, 'w') as out:
        out.write(f"{sample}\tCoding vs non coding\t")
        out.write(f"[coding bp, total bp, fraction coding]\t")
        out.write(f"{protlen}\t{dnalen}\t{protlen/dnalen}\n")





rule all:
    input:
        expand(os.path.join(STATS, "{sample}.average.phage.fraction"), sample=SAMPLES),
        os.path.join(STATS, "all_stats.tsv")


rule call_orfs:
    input:
        os.path.join(CONTIGS, "{sample}.fasta")
    output:
        os.path.join(ORFS, "{sample}.orfs.fna")
    shell:
        "phanotate.py -f fasta -o {output} {input}"

rule translate_orfs:
    """
    transeq is part of the emboss package
    Note: the -trim option removes trailing * or X at the end of the
    sequence
    """

    input:
        os.path.join(ORFS, "{sample}.orfs.fna")
    output:
        os.path.join(ORFS, "{sample}.orfs.faa")
    shell:
        "transeq -trim {input} {output}"

rule blast_phage_proteins:
    """
    This requires a phages.faa fasta database
    that we do not yet provide
    """
    input:
        os.path.join(ORFS, "{sample}.orfs.faa")
    output:
        os.path.join(BLAST, "{sample}.phages.blastp")
    threads:
        8
    shell:
        "blastp -query {input} -db {PHAGEDB} -outfmt '6 std qlen slen' -out {output} -num_threads {threads}"

rule blast_nr:
    """
    This requires an nr database that you can download from NCBI
    """
    input:
        os.path.join(ORFS, "{sample}.orfs.faa")
    output:
        os.path.join(BLAST, "{sample}.nr.blastp")
    params:
        db = config['databases']['nr']
    threads:
        8
    shell:
        "blastp -query {input} -db {params.db} -outfmt '6 std qlen slen' -out {output} -num_threads {threads}"


rule average_phage_protein_len:
    input:
        bp = os.path.join(BLAST, "{sample}.phages.blastp")
    output:
        fr = os.path.join(BLAST, "{sample}.phage_prot_fractions.tsv"),
        st = os.path.join(STATS, "{sample}.average.phage.fraction")
    params:
        sample = "{sample}"
    run:
        av_protein_lengths(params.sample, input.bp, output.fr, output.st, 'phage')

rule average_nr_protein_len:
    input:
        bp = os.path.join(BLAST, "{sample}.nr.blastp")
    output:
        fr = os.path.join(BLAST, "{sample}.nr_prot_fractions.tsv"),
        st = os.path.join(STATS, "{sample}.average.nr.fraction")
    params:
        sample = "{sample}"
    run:
        av_protein_lengths(params.sample, input.bp, output.fr, output.st, 'nr')


rule adjacent_phage_orfs_same_protein:
    input:
        fa = os.path.join(ORFS, "{sample}.orfs.faa"),
        bp = os.path.join(BLAST, "{sample}.phages.blastp")
    output:
        st = os.path.join(STATS, "{sample}.phage.adjacent.tsv"),
        nh = os.path.join(STATS, "{sample}.phage.nohits.tsv")
    params:
        sample = "{sample}"
    run:
        count_adjacent_orfs(params.sample, input.fa, input.bp, output.st, output.nh, "phage")


rule adjacent_nr_orfs_same_protein:
    input:
        fa = os.path.join(ORFS, "{sample}.orfs.faa"),
        bp = os.path.join(BLAST, "{sample}.nr.blastp")
    output:
        st = os.path.join(STATS, "{sample}.nr.adjacent.tsv"),
        nh = os.path.join(STATS, "{sample}.nr.nohits.tsv")
    params:
        sample = "{sample}"
    run:
        count_adjacent_orfs(params.sample, input.fa, input.bp, output.st, output.nh, "nr")


rule coding_vs_noncoding:
    input:
        ct = os.path.join(CONTIGS, "{sample}.fasta"),
        fa = os.path.join(ORFS, "{sample}.orfs.faa")
    output:
        cn = os.path.join(STATS, "{sample}.coding_noncoding.tsv")
    params:
        sample = "{sample}"
    run:
        coding_versus_noncoding(params.sample, input.ct, input.fa, output.cn)


rule combine_outputs:
    input:
        expand(os.path.join(STATS, "{sample}.average.phage.fraction"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.average.nr.fraction"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.phage.adjacent.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.phage.nohits.tsv"), sample=SAMPLES),
        #expand(os.path.join(STATS, "{sample}.nr.adjacent.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.coding_noncoding.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.nr.nohits.tsv"), sample=SAMPLES)
    output:
        os.path.join(STATS, "all_stats.tsv")
    shell:
        "cat {input}  > {output}"
