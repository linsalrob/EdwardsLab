"""
Snakefile to upload all the genomes to PATRIC
"""


# where is the data
FASTADIR = config['fasta']
OUTPUTDIR = config['output']


FASTA, = glob_wildcards(os.path.join(FASTADIR, '{fasta}.fasta'))

rule all:
    input:
        

