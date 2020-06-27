

import os
import sys
import subprocess


phispydir = "/home3/redwards/GitHubs/PhiSpy/test_genbank_files"
GENOMES, = glob_wildcards(os.path.join(phispydir, '{genome}.gb.gz'))


# the testing sets are "data/trainSet_{ID}.txt"
TRAININGS = ["32002.17", "272558.23", "224308.360", "411479.31", "206672.37", "224914.79", "190650.21", "195102.53", "212717.31", "243230.96", "1351.557", "199310.168", "83333.998", "83334.295", "155864.289", "71421.45", "272623.42", "272626.22", "169963.176", "266835.41", "83331.121", "83332.460", "122586.26", "122587.18", "272843.53", "208964.452", "160488.79", "267608.42", "220341.87", "211586.69", "198214", "1280.10152", "196620.15", "158878.38", "160490.61", "198466.10", "186103.26", "243277.252", "190486.46", "160492.65", "183190.38", "214092.200", "187410.24", "1367847.3", "318586.5", "1525717.3", "1660154.3", "147645.106"]


rule all:
    input:
        expand(os.path.join("test.v.train", "{genome}.phispy.train.{train}", "tptn.txt"),
               genome=GENOMES, train=TRAININGS)

rule run_phispy:
    input:
        g = os.path.join(phispydir, "{genome}.gb.gz")
    params:
        t = "data/trainSet_{train}.txt",
        o = os.path.join("test.v.train", "{genome}.phispy.train.{train}")
    output:
        temporary(os.path.join("test.v.train", "{genome}.phispy.train.{train}", "bacteria.fasta")),
        temporary(os.path.join("test.v.train", "{genome}.phispy.train.{train}", "bacteria.gbk")),
        temporary(os.path.join("test.v.train", "{genome}.phispy.train.{train}", "phage.fasta")),
        os.path.join("test.v.train", "{genome}.phispy.train.{train}", "phage.gbk"),
        os.path.join("test.v.train", "{genome}.phispy.train.{train}", "phispy.log"),
    shell:
        """
        PhiSpy.py -t {params.t}  -o {params.o} --output_choice 4 {input.g}
        """

rule count_tp_tn:
    input:
        gen = os.path.join(phispydir, "{genome}.gb.gz"),
        phg = os.path.join("test.v.train", "{genome}.phispy.train.{train}", "phage.gbk")
    output:
        tp = os.path.join("test.v.train", "{genome}.phispy.train.{train}", "tptn.txt")
    shell:
        """
        python3 ~/GitHubs/PhiSpy/scripts/compare_predictions_to_phages.py -t {input.gen} -p {input.phg} > {output.tp}
        """
