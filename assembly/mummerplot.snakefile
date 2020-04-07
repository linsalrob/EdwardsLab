"""
A snakefile to make a set of all vs. all pairwise
comparisons using mummer, and plot the images.

To run this snakefile you will need:

    - [mummer](http://mummer.sourceforge.net/)
    - [snakemake](https://snakemake.readthedocs.io/)
    - [imagemagick](https://imagemagick.org/index.php) for the montage command
    - probably gnuplot

To run the command, you will need the mummerplot.yaml. I recommend
you copy that to your local directory and use that version of it.

Then you can run this file with

snakemake --configfile mummerplot.yaml mummerplot.snakefile

"""

import os
import re

# please define your inputs and outputs in this
# yaml file. That way, you can copy that file
# to any directory, and run this same
# snakefile

# please then define it on the command line as above
# configfile: "mummerplot.yaml"

## User configurable options
# directory with the fasta files, one per genome
# currently we do not support gzipped fasta files
# though I could add that if required.
FASTADIR = config['fasta']
# this is a directory for all the reverse complemented fasta files (as necessary)
FASTARCDIR = config['fasta_rc']
# ouput directory for the mummer intermediary files
OUTDIR   = config['mummer_output']
# output directory for the images from mummerplot.
PNGDIR = config["mummer_png"]
# the filename for the final image that will be a montage
# of all the mummer plots

outputfilename = config['montage']

# The only bit that is tricky is that your fasta filenames
# - must end .fasta (not .fna or .fa)
# - must (currently) not have any periods in them.
# There is an issue with how snakemake handles
# wildcards that I am still resolving, and
# this is the simplest solution!

REVCOM = os.path.join(os.environ["HOME"], "GitHubs/EdwardsLab/mummer/reverse_complement_fasta.py")


# define fastafilename to be only ^\w+.fasta$. This is important to stop cyclic 
# use of the wildcard as files are created.
wildcard_constraints:
    faf1 = '\w+',
    faf2 = '\w+'

# you will need this in your path somewhere


# mummer needs to know the length of the genome
# to appropriately set the x- and y- axis of the plot
# figure out the fasta files and the size of the genomes

genome_len = {}
filename = {}
def fasta_len(faf):
    """
    calculate the length of sequence in the file
        - we pass this to mummer
        - we also filter to remove any empty fasta files 
               (yes, my assembly pipeline can make those!)
    Note that in this step, we trim the .fasta, .fna, or .fa
    extension off of the file, and so we also store the file
    name to use in the first step
    """
    seqlen = 0
    global genome_len
    global filename
    with open(os.path.join(FASTADIR, f"{faf}.fasta"), 'r') as f:
        for l in f:
            if not l.startswith('>'):
                seqlen += len(l.strip())
    bc = faf.replace('.fasta', '')
    filename[bc] = faf
    genome_len[bc] = seqlen
    print(f"GENOME LEN {faf} {bc} {seqlen}")
    return seqlen


# a simple function to convert a wildcard into a filename
# we just have to do this to dereference the hash
def wc2fn1(wildcards):
    return filename[wildcards.faf1]
def wc2fn2(wildcards):
    return filename[wildcards.faf2]

# find the fasta files
FASTA, = glob_wildcards(os.path.join(FASTADIR, '{fasta}.fasta'))
# if your files end .fna use this line instead
# FASTA, = glob_wildcards(os.path.join(FASTADIR, '{fasta}.fna'))

# now calculate the lengths of those sequences
lengths = list(map(fasta_len, FASTA))

# now find the non empty fasta files and remove them
FA = [fa[0] for fa in filter(lambda g: g[1] > 0, genome_len.items())]
FA.sort()


# this is the rule that calls everything
rule all:
    input:
        # this figures out the rc as necessary
        expand(os.path.join(FASTARCDIR, "{sample}.fasta"), sample=FA),
        # this runs mummer and then mummerplot of all fasta files against each other
        expand(os.path.join(PNGDIR, "{faf1}.{faf2}.png"), faf1 = FA, faf2 = FA),
        # this creates the final montage
        os.path.join(PNGDIR, outputfilename)



rule reverse_complement:
    """
    test if we need to reverse complement some sequences
    """
    input:
        expand(os.path.join(FASTADIR, "{sample}.fasta"), sample=FA)
    output:
        expand(os.path.join(FASTARCDIR, "{sample}.fasta"), sample=FA)
    shell:
        """
        python3 {REVCOM} -d {FASTADIR} -k 6 -o {FASTARCDIR} -v
        """


rule run_mummer:
    """
    Here, we run mummer on each of the files.
    Remember, we trimmed the extension off earlier, so here we need the
    real file name.
    """
    input:
        ref = os.path.join(FASTARCDIR, "{faf1}.fasta"),
        test = os.path.join(FASTARCDIR, "{faf2}.fasta")
    output:
        mums = os.path.join(OUTDIR, "{faf1}.{faf2}.mums")
    shell:
        "mummer -mum -b -c {input.ref} {input.test} > {output.mums}"


rule make_plot:
    """
    Take the mummer output and make a plot of it.
    gp is the gnuplot file
    png is the png file

    There is a little bit of magic in this, because
    if the file doesn't exist (i.e. the sequence 
    was empty or there are no mums between the two
    sequences), we just make a white square
    for that space. This happens when mummerplot
    exits with a code of 25.

    We also add a label to the png file that
    is used in the montage.
    See https://www.imagemagick.org/Usage/montage/#label
    """
    input:
        mums = os.path.join(OUTDIR, "{faf1}.{faf2}.mums")
    output:
        fp  = os.path.join(OUTDIR, "{faf1}.{faf2}.gp"),
        png = os.path.join(PNGDIR, "{faf1}.{faf2}.png")
    params:
        g1 = lambda wildcards: genome_len[wildcards.faf1],
        g2 = lambda wildcards: genome_len[wildcards.faf2],
        label = "{faf1} {faf2}",
        outbase = os.path.join(OUTDIR, "{faf1}.{faf2}")
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
    """
    Combine all the mummerplots into a single montage
    """
    input:
        sorted(expand(os.path.join(PNGDIR, "{bc1}.{bc2}.png"), bc1 = FA, bc2 = FA))
    output:
        os.path.join(PNGDIR, outputfilename)
    shell:
        "montage {input} {output}"
