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
        i.      Number of contigs
        ii.     Length
        iii.    %GC
        iv.     stdev of %GC in 20bp windows
        v.      Fraction of “N” bases
        vi.     N50 (add L50/U50?)
        vii.    N75
        viii.   How many frameshifted proteins
        ix.     genome length too large
        x.      genome length too small
        xi.     abnormal gene to sequence ratio
        xii.    low gene count
        xiii.   average gene length
        ixx.    compare using translated blast (blastp or diamond??) and check for same hits on
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
import socket
from roblib import bcolors, stream_fasta, mean, stdev, median, is_hypothetical
from pppf_databases import connect_to_db, disconnect
from pppf_clusters import proteinid_to_function

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()

sys.stderr.write(f"Running on {socket.gethostname()}\n")

CONTIGS = config['output_paths']['contigs']
ORFS    = config['output_paths']['orfs']
BLAST   = config['output_paths']['blast']
STATS   = config['output_paths']['statistics']
RESULTS = config['output_paths']['results']

# executables
BLASTEXC = config['executable_paths']['blast']

# our phage and bacterial protein databases
PHAGEDB = os.path.join(config['install_paths']['databases'], config['databases']['phage_proteins'])
if not os.path.exists(PHAGEDB):
    sys.stderr.write(f"FATAL: {PHAGEDB} not found. Please check your paths\n")
    sys.exit()
BACTDB = os.path.join(config['install_paths']['databases'], config['databases']['bacterial_proteins'])
if not os.path.exists(BACTDB):
    sys.stderr.write(f"FATAL: {BACTDB} not found. Please check your paths\n")
    sys.exit()

# include functions from phage clusters
PHAGECLUSTERS = os.path.join(config['databases']['phage_cluster_database'])
if not os.path.exists(PHAGECLUSTERS):
    sys.stderr.write(f"Phage Cluster database {PHAGECLUSTERS} not found. Skipping functional analysis\n")

SAMPLES, = glob_wildcards(os.path.join(CONTIGS, '{sample}.fasta'))


def kmer_kurtosis(sequences, k=7):
    """
    Calculate the kurtosis of k-mers in the sequences
    :param sequences: dict of sequences
    :param k: kmer size to use
    """

    # count the abundance of all kmers
    kmers = {}
    for s in sequences:
        i=0
        while i < len(sequences[s])-k:
            subseq = sequences[s][i:i+k]
            kmers[subseq] = kmers.get(subseq, 0) + 1
            i += 1
    counts = list(kmers.values())
    av = mean(counts)
    st = stdev(counts)

    kurtosis = 0
    for i,j in enumerate(counts):
        kurtosis += (j-av)**4

    kurtosis = kurtosis / (len(counts) * (st**4))
    kurtosis -= 3

    return kurtosis



def contig_stats(sample, fastafile, contigstats):
    """
    Calculate some stats about the contigs in the fasta file
    """

    seqs = {}
    s = ""
    oldseqid = None
    with open(fastafile, 'r') as f:
        for l in f:
            if l.startswith('>'):
                p = l.split()
                seqid = p[0].replace('>', '', 1)
                if oldseqid:
                    seqs[oldseqid]=s
                    s = ""
                oldseqid = seqid
            else:
                s += l.strip()
    seqs[oldseqid]=s

    lens = [len(seqs[x]) for x in seqs]
    gc = 0
    fn = 0
    for s in seqs:
        gc += seqs[s].upper().count("G")
        gc += seqs[s].upper().count("C")
        fn += seqs[s].upper().count("N")

    gc /= sum(lens)
    fn /= sum(lens)


    with open(contigstats, "w") as out:
        out.write(f"{sample}\tNumber of contigs\t{len(seqs)}\n")
        out.write(f"{sample}\tTotal length\t{sum(lens)}\n")
        out.write(f"{sample}\tLongest sequence\t{max(lens)}\n")
        out.write(f"{sample}\tk-mer kurtosis\t{kmer_kurtosis(seqs, 7)}\n")
        out.write(f"{sample}\tpercent GC\t{gc}\n")
        out.write(f"{sample}\tFraction of N's\t{fn}")




def count_adjacent_orfs(sample, fastafile, blastfile, adjacentout, nohitsout, searchtype):
    """
    In this version, we just count how many query proteins map to the
    same target protein. This is easier than the original approach 
    below
    """
    sys.stderr.write(f"{bcolors.GREEN}Counting adjacent ORFs for {sample} and {searchtype}{bcolors.ENDC}\n")
    orfs = set()
    with open(fastafile, 'r') as f:
        for l in f:
            if l.startswith('>'):
                p = l.split()
                seqid = p[0].replace('>', '', 1)
                orfs.add(seqid)
    
    hits = {}
    queries = set()
    with open(blastfile, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if p[1] not in hits:
                hits[p[1]] = set()
            hits[p[1]].add(p[0])
            queries.add(p[0])
    count = 0
    for h in hits:
        if len(hits[h]) > 1:
            count += 1
    
    nohits = (len(orfs) - len(queries))
    
    
    with open(adjacentout, 'w') as out:
        out.write(f"{sample}\tNumber of targets with two {searchtype} queries\t")
        out.write(f"[two queries, total orfs, fraction twos]\t")
        out.write(f"{count}\t{len(orfs)}\t{count/len(orfs)}\n")
    with open(nohitsout, 'w') as out:
        out.write(f"{sample}\tNumber of orfs with no hits\t")
        out.write(f"[nohits, total orfs, fraction no hits]\t")
        out.write(f"{nohits}\t{len(orfs)}\t{nohits/len(orfs)}\n")


def count_adjacent_orfs_original(sample, fastafile, blastfile, adjacentout, nohitsout, searchtype):
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
        out.write("[num orfs, median proportional length, average proportional length]\t")
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

def check_phage_functions(sample, blastfile, outputfile):
    """
    Count how many proteins have hypothetical functions
    """

    out = open(outputfile, 'w')
    if not os.path.exists(PHAGECLUSTERS):
        out.close()
        return
    
    phage_cluster_db = connect_to_db(PHAGECLUSTERS)
    hypo = 0
    nonhypo = 0
    with open(blastfile, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            fn = proteinid_to_function(p[1], phage_cluster_db)
            if is_hypothetical(fn):
                hypo +=1
            else:
                nonhypo += 1
    out.write(f"{sample}\tHypothetical proteins\t")
    out.write("[Hypothetical, Non-hypothetical, Fraction hypothetical]\t")
    out.write(f"{hypo}\t{nonhypo}\t{hypo/(hypo+nonhypo)}\n")
    out.close()

        


rule all:
    input:
        expand(os.path.join(STATS, "{sample}.average.phage.fraction.tsv"), sample=SAMPLES),
        os.path.join(RESULTS, "all_stats.tsv"),
        os.path.join(RESULTS, "average.phage.fraction.tsv"),
        os.path.join(RESULTS, "phage.hyponon.tsv"),
        os.path.join(RESULTS, "bacteria.nohits.tsv"),
        os.path.join(RESULTS, "coding_noncoding.tsv"),
        os.path.join(RESULTS, "bacteria.adjacent.tsv"),
        os.path.join(RESULTS, "phage.nohits.tsv"),
        os.path.join(RESULTS, "phage.adjacent.tsv"),
        os.path.join(RESULTS, "average.bacteria.fraction.tsv")


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
    params:
        blastp = os.path.join(BLASTEXC, "blastp")
    threads:
        8
    shell:
        "{params.blastp} -query {input} -db {PHAGEDB} -outfmt '6 std qlen slen' -out {output} -num_threads {threads}"

rule blast_bact_proteins:
    """
    This requires a phages.faa fasta database
    that we do not yet provide
    """
    input:
        os.path.join(ORFS, "{sample}.orfs.faa")
    output:
        os.path.join(BLAST, "{sample}.bacteria.blastp")
    params:
        blastp = os.path.join(BLASTEXC, "blastp")
    threads:
        8
    shell:
        "{params.blastp} -query {input} -db {BACTDB} -outfmt '6 std qlen slen' -out {output} -num_threads {threads}"

rule blast_nr:
    """
    This requires an nr database that you can download from NCBI
    """
    input:
        os.path.join(ORFS, "{sample}.orfs.faa")
    output:
        os.path.join(BLAST, "{sample}.nr.blastp")
    params:
        blastp = os.path.join(BLASTEXC, "blastp"),
        db = config['databases']['nr']
    threads:
        8
    shell:
        "{params.blastp} -query {input} -db {params.db} -outfmt '6 std qlen slen' -out {output} -num_threads {threads}"

rule average_phage_protein_len:
    input:
        bp = os.path.join(BLAST, "{sample}.phages.blastp")
    output:
        fr = os.path.join(BLAST, "{sample}.phage_prot_fractions.tsv"),
        st = os.path.join(STATS, "{sample}.average.phage.fraction.tsv")
    params:
        sample = "{sample}",
    run:
        av_protein_lengths(params.sample, input.bp, output.fr, output.st, 'phage')

rule average_bact_protein_len:
    input:
        bp = os.path.join(BLAST, "{sample}.bacteria.blastp")
    output:
        fr = os.path.join(BLAST, "{sample}.bact_prot_fractions.tsv"),
        st = os.path.join(STATS, "{sample}.average.bacteria.fraction.tsv")
    params:
        sample = "{sample}"
    run:
        av_protein_lengths(params.sample, input.bp, output.fr, output.st, 'bacteria')

rule average_nr_protein_len:
    input:
        bp = os.path.join(BLAST, "{sample}.nr.blastp")
    output:
        fr = os.path.join(BLAST, "{sample}.nr_prot_fractions.tsv"),
        st = os.path.join(STATS, "{sample}.average.nr.fraction.tsv")
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

rule adjacent_bacteria_orfs_same_protein:
    input:
        fa = os.path.join(ORFS, "{sample}.orfs.faa"),
        bp = os.path.join(BLAST, "{sample}.bacteria.blastp")
    output:
        st = os.path.join(STATS, "{sample}.bacteria.adjacent.tsv"),
        nh = os.path.join(STATS, "{sample}.bacteria.nohits.tsv")
    params:
        sample = "{sample}"
    run:
        count_adjacent_orfs(params.sample, input.fa, input.bp, output.st, output.nh, "bacteria")


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

rule hypo_vs_non:
    input:
        bf = os.path.join(BLAST, "{sample}.phages.blastp")
    output:
        hn = os.path.join(STATS, "{sample}.phage.hyponon.tsv")
    params:
        sample = "{sample}"
    run:
        check_phage_functions(params.sample, input.bf, output.hn)

rule count_contig_stats:
    input:
        ct = os.path.join(CONTIGS, "{sample}.fasta")
    output:
        cs = os.path.join(STATS, "{sample}.contig_stats.tsv")
    params:
        sample = "{sample}"
    run:
        contig_stats(params.sample, input.ct, output.cs)

rule combine_average_phage_fraction:
    input:
        expand(os.path.join(STATS, "{sample}.average.phage.fraction.tsv"), sample=SAMPLES),
    output:
        os.path.join(RESULTS, "average.phage.fraction.tsv")
    shell:
        "cat {input} > {output}"

rule combine_phage_hyponon:
    input:
        expand(os.path.join(STATS, "{sample}.phage.hyponon.tsv"), sample=SAMPLES)
    output:
        os.path.join(RESULTS, "phage.hyponon.tsv")
    shell:
        "cat {input} > {output}"

rule combine_bacteria_nohits:
    input:
        expand(os.path.join(STATS, "{sample}.bacteria.nohits.tsv"), sample=SAMPLES),
    output:
        os.path.join(RESULTS, "bacteria.nohits.tsv")
    shell:
        "cat {input} > {output}"

rule combine_coding_noncoding:
    input:
        expand(os.path.join(STATS, "{sample}.coding_noncoding.tsv"), sample=SAMPLES),
    output:
        os.path.join(RESULTS, "coding_noncoding.tsv")
    shell:
        "cat {input} > {output}"

rule combine_bacteria_adjacent:
    input:
        expand(os.path.join(STATS, "{sample}.bacteria.adjacent.tsv"), sample=SAMPLES),
    output:
        os.path.join(RESULTS, "bacteria.adjacent.tsv")
    shell:
        "cat {input} > {output}"

rule combine_phage_nohits:
    input:
        expand(os.path.join(STATS, "{sample}.phage.nohits.tsv"), sample=SAMPLES),
    output:
        os.path.join(RESULTS, "phage.nohits.tsv")
    shell:
        "cat {input} > {output}"

rule combine_phage_adjacent:
    input:
        expand(os.path.join(STATS, "{sample}.phage.adjacent.tsv"), sample=SAMPLES),
    output:
        os.path.join(RESULTS, "phage.adjacent.tsv")
    shell:
        "cat {input} > {output}"

rule combine_average_bacteria_fraction:
    input:
        expand(os.path.join(STATS, "{sample}.average.bacteria.fraction.tsv"), sample=SAMPLES),
    output:
        os.path.join(RESULTS, "average.bacteria.fraction.tsv")
    shell:
        "cat {input} > {output}"

rule combine_contig_stats:
    input:
        expand(os.path.join(STATS, "{sample}.contig_stats.tsv"), sample=SAMPLES)
    output:
        os.path.join(RESULTS, "contig_stats.tsv")
    shell:
        "cat {input} > {output}"
        
rule combine_outputs:
    """
    This rule does two things:
        1. We make sure all the samples are run
        2. We combine all the results, but they are already combined above, hopefully!
    """
    input:
        os.path.join(RESULTS, "average.bacteria.fraction.tsv"),
        os.path.join(RESULTS, "average.phage.fraction.tsv"),
        os.path.join(RESULTS, "bacteria.adjacent.tsv"),
        os.path.join(RESULTS, "bacteria.nohits.tsv"),
        os.path.join(RESULTS, "coding_noncoding.tsv"),
        os.path.join(RESULTS, "contig_stats.tsv"),
        os.path.join(RESULTS, "phage.adjacent.tsv"),
        os.path.join(RESULTS, "phage.hyponon.tsv"),
        os.path.join(RESULTS, "phage.nohits.tsv"),

        
        # expand(os.path.join(STATS, "{sample}.average.phage.fraction.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.average.bacteria.fraction.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.bacteria.adjacent.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.bacteria.nohits.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.coding_noncoding.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.phage.adjacent.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.phage.nohits.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.phage.hyponon.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.contig_stats.tsv"), sample=SAMPLES),
        # uncomment these to run the nr blast, but it takes  a long time!
        # expand(os.path.join(STATS, "{sample}.average.nr.fraction"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.nr.adjacent.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.nr.nohits.tsv"), sample=SAMPLES),
    output:
        os.path.join(RESULTS, "all_stats.tsv")
    shell:
        "cat {input} | sort > {output}"
