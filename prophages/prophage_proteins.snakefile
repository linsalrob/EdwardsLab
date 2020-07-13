"""
A snakefile to download all genomes from RefSeq bacteria and run
PhiSpy on those genomes, create a phage genbank file, and then list
protein IDs that are phages

Note: before using this ensure rsync is not aliased!

"""

import os
import sys
from Bio import SeqIO
import gzip

GBKDIR = "genbank"
PHSDIR = "phispy"
PHGDIR = "phage_proteins"
STATSDIR = "stats"

# we use assembly_summary.tsv for the list of genomes, which is available
# from  ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
# however, we set that as a config so you can override it

assf = "assembly_summary.tsv"
if 'assemblies' in config:
    sys.stderr.write(f"Using {config['assemblies']} as the assembly_summary file\n")
    assf = config['assemblies']


def get_samples():
    samples = {}
    with open(assf, 'r') as f:
        for l in f:
            if l.startswith("#"):
                continue
            p = l.strip().split("\t")
            # we strip off the protocol because we prefer rsync but fall
            # fall back to curl if that doesn't work
            if p[19] == "na":
                continue
            p[19] = p[19].replace("ftp://", "", 1)
            ass_id = p[19].split("/")[-1]
            samples[ass_id] = p[19]

    return list(samples.keys()), samples

SAMPLES, URLS = get_samples()

def get_url(wildcards):
    return f"{URLS[wildcards.sample]}/{wildcards.sample}_genomic.gbff.gz"


def contig_lens(inf, outf):
    """
    Just report the contig lengths
    """
    with open(outf, 'w') as out:
        handle = gzip.open(inf, 'rt')
        for seq in SeqIO.parse(handle, "genbank"):
            out.write(f"{seq.id}\t{len(seq.seq)}\n")

rule all:
    input:
        expand(os.path.join(PHSDIR, "{sample}_prophage_coordinates.tsv"), sample=SAMPLES),
        expand(os.path.join(STATSDIR, "{sample}_contigs.txt"), sample=SAMPLES),
        expand(os.path.join(PHGDIR, "{sample}_protein_functions.txt"), sample=SAMPLES)

rule download_genbank:
    """
    We start with rsync (it is the preferred protocol per 
    https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#protocols
    but that seems to fail quite frequently (perhaps gets swamped by
    multiple requests?) and so we fall back to curl/https.

    Note that Rob has an alias to rsync to use ssh by default, 
    so unaliasing that unsets it
    """
    output:
        #temporary(os.path.join(GBKDIR, "{sample}_genomic.gbff.gz"))
        os.path.join(GBKDIR, "{sample}_genomic.gbff.gz")
    params:
        url = get_url,
        gbd = GBKDIR
    shell:
        """
        set +e
        unalias rsync
        rsync --copy-links --recursive --times --verbose rsync://{params.url} {params.gbd}
        exitcode=$?
        if [ $exitcode == 10 ]
        then
            curl -Lo {output} https://{params.url}
        fi
        """
       
"""
        bfa = temporary(os.path.join(PHSDIR, "{sample}_bacteria.fasta")),
        pfa = temporary(os.path.join(PHSDIR, "{sample}_phage.fasta")),
        lgf = temporary(os.path.join(PHSDIR, "{sample}_phispy.log")),
        bct = temporary(os.path.join(PHSDIR, "{sample}_bacteria.gbk")),
        gbk = temporary(os.path.join(PHSDIR, "{sample}_phage.gbk")),
        tsv = os.path.join(PHSDIR, "{sample}_prophage_coordinates.tsv"),
"""

rule run_phispy:
    input:
        os.path.join(GBKDIR, "{sample}_genomic.gbff.gz")
    output:
        temporary(os.path.join(PHSDIR, "{sample}_bacteria.fasta")),
        temporary(os.path.join(PHSDIR, "{sample}_phage.fasta")),
        os.path.join(PHSDIR, "{sample}_phispy.log"),
        temporary(os.path.join(PHSDIR, "{sample}_bacteria.gbk")),
        os.path.join(PHSDIR, "{sample}_phage.gbk"),
        os.path.join(PHSDIR, "{sample}_prophage_coordinates.tsv"),
    params:
        phispydir = PHSDIR,
        sample = "{sample}"
    benchmark:
        "benchmarks/{sample}.phispy.txt"
    shell:
        """
        set +e
        PhiSpy.py -o {params.phispydir} --quiet --output_choice 5 -p {params.sample} {input}
        exitcode=$?
        if [ $exitcode -gt 19 ]
        then
            for o in {output}; do touch $o; done
        fi
        """

rule create_phage_protein_list:
    input:
        gbk = os.path.join(PHSDIR, "{sample}_phage.gbk")
    output:
        txt = os.path.join(PHGDIR, "{sample}_protein_functions.txt")
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/bin/genbank2sequences.py -g {input.gbk} --functions {output.txt}
        """
    

rule contig_stats:
    input:
        gbk = os.path.join(GBKDIR, "{sample}_genomic.gbff.gz")
    output:
        out = os.path.join(STATSDIR, "{sample}_contigs.txt")
    run:
        contig_lens(input.gbk, output.out)


