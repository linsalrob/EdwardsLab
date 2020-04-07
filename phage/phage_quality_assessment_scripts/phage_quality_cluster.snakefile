"""
A pipeline to check the quality of phage contigs.

Note: In this version I refactored all the run: sections to scripts so that they get run independently on the cluster
We were having a delay waiting for the run's to finish!
Note that we leave the phage functions because that uses SQL so we can't parallelize it.

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
from roblib import bcolors, stream_fasta, median, is_hypothetical
from pppf_databases import connect_to_db, disconnect
from pppf_clusters import proteinid_to_function

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()

CONTIGS = config['paths']['contigs']
ORFS    = config['paths']['orfs']
BLAST   = config['paths']['blast']
STATS   = config['paths']['statistics']

# our phage and bacterial protein databases
PHAGEDB = os.path.join(config['paths']['databases'], config['databases']['phage_proteins'])
if not os.path.exists(PHAGEDB):
    sys.stderr.write(f"FATAL: {PHAGEDB} not found. Please check your paths\n")
    sys.exit()
BACTDB = os.path.join(config['paths']['databases'], config['databases']['bacterial_proteins'])
if not os.path.exists(BACTDB):
    sys.stderr.write(f"FATAL: {BACTDB} not found. Please check your paths\n")
    sys.exit()

# include functions from phage clusters
PHAGECLUSTERS = os.path.join(config['databases']['phage_cluster_database'])
if not os.path.exists(PHAGECLUSTERS):
    sys.stderr.write(f"Phage Cluster database {PHAGECLUSTERS} not found. Skipping functional analysis\n")

SAMPLES, = glob_wildcards(os.path.join(CONTIGS, '{sample}.fasta'))




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
    disconnect(phage_cluster_db)




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

rule blast_bact_proteins:
    """
    This requires a phages.faa fasta database
    that we do not yet provide
    """
    input:
        os.path.join(ORFS, "{sample}.orfs.faa")
    output:
        os.path.join(BLAST, "{sample}.bacteria.blastp")
    threads:
        8
    shell:
        "blastp -query {input} -db {BACTDB} -outfmt '6 std qlen slen' -out {output} -num_threads {threads}"

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
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/phage/phage_quality_assessment_scripts/av_protein_lengths.py \
        -s {params.sample} -b {input.bp} -f {output.fr} -o {output.st} -t "phage"
        """

rule average_bact_protein_len:
    input:
        bp = os.path.join(BLAST, "{sample}.bacteria.blastp")
    output:
        fr = os.path.join(BLAST, "{sample}.bact_prot_fractions.tsv"),
        st = os.path.join(STATS, "{sample}.average.bacteria.fraction")
    params:
        sample = "{sample}"
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/phage/phage_quality_assessment_scripts/av_protein_lengths.py \
        -s {params.sample} -b {input.bp} -f {output.fr} -o {output.st} -t "bacteria"
        """

rule average_nr_protein_len:
    input:
        bp = os.path.join(BLAST, "{sample}.nr.blastp")
    output:
        fr = os.path.join(BLAST, "{sample}.nr_prot_fractions.tsv"),
        st = os.path.join(STATS, "{sample}.average.nr.fraction")
    params:
        sample = "{sample}"
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/phage/phage_quality_assessment_scripts/av_protein_lengths.py \
        -s {params.sample} -b {input.bp} -f {output.fr} -o {output.st} -t "nr"
        """


rule adjacent_phage_orfs_same_protein:
    input:
        fa = os.path.join(ORFS, "{sample}.orfs.faa"),
        bp = os.path.join(BLAST, "{sample}.phages.blastp")
    output:
        st = os.path.join(STATS, "{sample}.phage.adjacent.tsv"),
        nh = os.path.join(STATS, "{sample}.phage.nohits.tsv")
    params:
        sample = "{sample}"
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/phage/phage_quality_assessment_scripts/count_adjacent_orfs.py \
        -s {params.sample} -f {input.fa} -b {input.bp} -a {output.st} -n {output.nh} -t "phage"
        """

rule adjacent_bacteria_orfs_same_protein:
    input:
        fa = os.path.join(ORFS, "{sample}.orfs.faa"),
        bp = os.path.join(BLAST, "{sample}.bacteria.blastp")
    output:
        st = os.path.join(STATS, "{sample}.bacteria.adjacent.tsv"),
        nh = os.path.join(STATS, "{sample}.bacteria.nohits.tsv")
    params:
        sample = "{sample}"
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/phage/phage_quality_assessment_scripts/count_adjacent_orfs.py \
        -s {params.sample} -f {input.fa} -b {input.bp} -a {output.st} -n {output.nh} -t "bacteria"
        """

rule adjacent_nr_orfs_same_protein:
    input:
        fa = os.path.join(ORFS, "{sample}.orfs.faa"),
        bp = os.path.join(BLAST, "{sample}.nr.blastp")
    output:
        st = os.path.join(STATS, "{sample}.nr.adjacent.tsv"),
        nh = os.path.join(STATS, "{sample}.nr.nohits.tsv")
    params:
        sample = "{sample}"
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/phage/phage_quality_assessment_scripts/count_adjacent_orfs.py \
        -s {params.sample} -f {input.fa} -b {input.bp} -a {output.st} -n {output.nh} -t "nr"
        """


rule coding_vs_noncoding:
    input:
        ct = os.path.join(CONTIGS, "{sample}.fasta"),
        fa = os.path.join(ORFS, "{sample}.orfs.faa")
    output:
        cn = os.path.join(STATS, "{sample}.coding_noncoding.tsv")
    params:
        sample = "{sample}"
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/phage/phage_quality_assessment_scripts/coding_vs_noncoding.py \
        -s {params.sample} -c {input.ct} -f {input.fa} -o {output.cn}
        """

rule hypo_vs_non:
    input:
        bf = os.path.join(BLAST, "{sample}.phages.blastp")
    output:
        hn = os.path.join(STATS, "{sample}.phage.hyponon.tsv")
    params:
        sample = "{sample}"
    run:
        check_phage_functions(params.sample, input.bf, output.hn)


rule combine_outputs:
    input:
        expand(os.path.join(STATS, "{sample}.average.phage.fraction"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.average.bacteria.fraction"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.phage.adjacent.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.phage.nohits.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.bacteria.adjacent.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.coding_noncoding.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.bacteria.nohits.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS, "{sample}.phage.hyponon.tsv"), sample=SAMPLES)
        # uncomment these to run the nr blast, but it takes  a long time!
        # expand(os.path.join(STATS, "{sample}.average.nr.fraction"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.nr.adjacent.tsv"), sample=SAMPLES),
        # expand(os.path.join(STATS, "{sample}.nr.nohits.tsv"), sample=SAMPLES),
    output:
        os.path.join(STATS, "all_stats.tsv")
    shell:
        "cat {input}  > {output}"
