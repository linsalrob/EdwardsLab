

import os
import sys
import subprocess


phispydir = "/home3/redwards/GitHubs/PhiSpy/test_genbank_files"
GENOMES, = glob_wildcards(os.path.join(phispydir, '{genome}.gb.gz'))

rule all:
    input:
        expand(os.path.join("{genome}.phispy.pg{phage_genes}", "tptn.txt"),
               genome=GENOMES, phage_genes=[0,1,2,3,4,5])

rule run_phispy:
    input:
        g = os.path.join(phispydir, "{genome}.gb.gz")
    params:
        pg = "{phage_genes}",
        o = "{genome}.phispy.pg{phage_genes}"
    output:
        temporary(os.path.join("{genome}.phispy.pg{phage_genes}", "bacteria.fasta")),
        temporary(os.path.join("{genome}.phispy.pg{phage_genes}", "bacteria.gbk")),
        temporary(os.path.join("{genome}.phispy.pg{phage_genes}", "phage.fasta")),
        os.path.join("{genome}.phispy.pg{phage_genes}", "phage.gbk"),
        os.path.join("{genome}.phispy.pg{phage_genes}", "phispy.log"),
    shell:
        """
        PhiSpy.py --phage_genes {params.pg}  -o {params.o} --output_choice 4 {input.g}
        """

rule count_tp_tn:
    input:
        gen = os.path.join(phispydir, "{genome}.gb.gz"),
        phg = os.path.join("{genome}.phispy.pg{phage_genes}", "phage.gbk")
    output:
        tp = os.path.join("{genome}.phispy.pg{phage_genes}", "tptn.txt")
    shell:
        """
        python3 ~/GitHubs/PhiSpy/scripts/compare_predictions_to_phages.py -t {input.gen} -p {input.phg} > {output.tp}
        """
