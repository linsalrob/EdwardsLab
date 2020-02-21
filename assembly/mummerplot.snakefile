"""
Snake file to make mummerplots for all the fasta files. 

We need to do a couple of things to read the fasta files and get their
lengths to make the mummer plots.
"""

import os
import re

# figure out the fasta files and the size of the genomes

genome_len = {}

ASSDIR = "assembly"
OUTDIR = "mummer"
PNGDIR = "mummer_png"

fastare = re.compile(r'\.fasta$')

# define barcode to be only ^\w+$. This is important to stop cyclic 
# use of the wildcard as files are created.
wildcard_constraints:
   barcode1 = '\w+',
   barcode2 = '\w+'


def fasta_len(faf):
    """ calculate the length of sequence in the file to ensure > 0"""
    seqlen = 0
    global genome_len
    with open(os.path.join(ASSDIR, faf), 'r') as f:
        for l in f:
            if not l.startswith('>'):
                seqlen += len(l.strip())
    print(f"GENOME LEN {faf} {seqlen}")
    bc = faf.replace('.fasta', '')
    genome_len[bc] = seqlen
    return seqlen


# find all the fasta files
def fasta_filter(faf):
    if not fastare.search(faf):
        return False
    return True


ALLFAFILES = list(filter(fasta_filter, os.listdir(ASSDIR)))
lengths = list(map(fasta_len, ALLFAFILES))

# now find the non empty fasta files
FA = [fa[0] for fa in filter(lambda g: g[1] > 0, genome_len.items())]

print(f"BC: {FA}")

rule all:
    input:
        expand(os.path.join(PNGDIR, "{barcode1}.{barcode2}.png"), barcode1 = FA, barcode2 = FA),
        os.path.join(PNGDIR, "montage.png")


rule run_mummer:
    input:
        ref = os.path.join(ASSDIR, "{barcode1}.fasta"),
        test = os.path.join(ASSDIR, "{barcode2}.fasta")
    output:
        mums = os.path.join(OUTDIR, "{barcode1}.{barcode2}.mums")
    shell:
        "mummer -mum -b -c {input.ref} {input.test} > {output.mums}"


rule make_plot:
    input:
        mums = os.path.join(OUTDIR, "{barcode1}.{barcode2}.mums")
    output:
        fp  = os.path.join(OUTDIR, "{barcode1}.{barcode2}.gp"),
        png = os.path.join(PNGDIR, "{barcode1}.{barcode2}.png")
    params:
        g1 = lambda wildcards: genome_len[wildcards.barcode1],
        g2 = lambda wildcards: genome_len[wildcards.barcode2],
        label = "{barcode1} {barcode2}",
        outbase = os.path.join(OUTDIR, "{barcode1}.{barcode2}")
    shell:
        """
        set +e
        mummerplot -x '[0,{params.g1}]' -y '[0,{params.g2}]' --png -p {params.outbase} {input.mums}
        exitcode=$?
        if [ $exitcode == 25 ]
        then
            touch {output.fp}
            convert -label {params.label} -size 100x100 xc:white {output.png}
            exit 0
        else
            convert -label {params.label} {params.outbase}.png {output.png}
            exit $exitcode
        fi
        """

rule montage:
    input:
        sorted(expand(os.path.join(PNGDIR, "{bc1}.{bc2}.png"), bc1 = FA, bc2 = FA))
    output:
        os.path.join(PNGDIR, "montage.png")
    shell:
        "montage {input} {output}"
