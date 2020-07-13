


import os

gbk = "phage_100_genbank"
faa = "proteins"
fna = "nucleotides"
PHAGES, = glob_wildcards(os.path.join(gbk, '{phage}_phage.gbk'))


rule all:
    input:
        expand(os.path.join(fna, '{phage}.fna'), phage=PHAGES)

rule gbk2faa:
    input:
        os.path.join(gbk, '{phage}_phage.gbk')
    output:
        faa = os.path.join(faa, '{phage}.faa'),
        fna = os.path.join(fna, '{phage}.fna')
    shell:
        """
        python3 ~/GitHubs/EdwardsLab/bin/genbank2sequences.py -g {input} -a {output.faa} -n {output.fna} -c
        """

