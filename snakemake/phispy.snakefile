"""
Run phispy on a bunch of genomes from patric
"""

import os
import sys
from random import shuffle

configfile: "phispy.yaml"

IDS=[]
with open("genome_ids.txt", 'r') as f:
    for l in f:
        IDS.append(l.strip())


# choose some ids at random
shuffle(IDS)
IDS = IDS[0:config['number_of_genomes']]



GBFDIR = config["directories"]["genbank_files"]
PHIDIR   = config["directories"]["phispy_files"]
GTODIR   = config["directories"]["gto_files"]


rule all:
    input:
        #expand(os.path.join(GBFDIR, "{sample}.gbf"), sample=IDS)
        expand(os.path.join(PHIDIR, "{sample}", "prophage.tbl"), sample=IDS)

rule download_gto:
    output:
        os.path.join(GTODIR, "{sample}.gto")
    params:
        genome_id = "{sample}"
    shell:
        "p3-gto {params.genome_id}; mv {params.genome_id}.gto {output}"


rule export_genbank:
    input:
        os.path.join(GTODIR, "{sample}.gto")
    output:
        os.path.join(GBFDIR, "{sample}.gbf")
    shell:
        "rast-export-genome -i {input} -o {output} genbank"


rule run_phispy:
    input:
        os.path.join(GBFDIR, "{sample}.gbf")
    output:
        os.path.join(PHIDIR, "{sample}", "prophage_coordinates.tsv"),
        os.path.join(PHIDIR, "{sample}", "prophage.gff3"),
        os.path.join(PHIDIR, "{sample}", "prophage.tbl"),
        os.path.join(PHIDIR, "{sample}", "prophage_tbl.tsv"),
    params:
        outdir = os.path.join(PHIDIR, "{sample}")
    shell:
        "PhiSpy.py -o {params.outdir} {input}"


