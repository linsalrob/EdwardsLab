# Assembly and comparisons

These are some snakefiles for assembling nanopore data and plotting pairwise comparisons using mummer.

To run any of these you will need snakemake installed. Typically you'll want to do this in conda or something, but
I just use `pip3 install --user snakemake`

# Assembly pipeline

To build the assembly you will need:
- minimap2
- miniasm
- filtlong
- minipolish

Plus a few of the `count` scripts in [../bin](../bin) like `countfasta.py`, `countfastq.py` and `countgfa.py`

*Note:* This is a completely untrusted and untested pipeline that I have built and you should almost certainly not use it!


To run the assembly script on my cluster (ha, you can't), I use this command:

```
snakemake  --snakefile ~/EdwardsLab/assembly/nanopore_phage_assembly.snakefile \
	--cluster 'qsub -q default -cwd -o sge_out -e sge_err -V' \
	-j 200 --latency-wait 60 \
	> snakemake.out 2>&1
```

This pipeline is supposed to be better for phage genomes sequenced using the Nanopore: see [this paper](https://f1000research.com/articles/8-2138/v1) by Ryan Wick and Kat Holt.

# Pairwise Mummer

To make a pairwise plot of all the genome assemblies against each other you will need to install

- [mummer](http://mummer.sourceforge.net/)
- [snakemake](https://snakemake.readthedocs.io/)
- [imagemagick](https://imagemagick.org/index.php) for the montage command
- probably gnuplot

Then your files:
- need to be in fasta format
- should be in a directory called `fasta`
- should be in files that end with `.fasta`. Sorry, I can't yet do `.fna` or `.fa` files, though you could easily edit this snakefile to make that happen
- should only have letters, numbers, and underscores (ie. [A-Za-z0-9\_] in the filenames). There is an issue with recursion that I also don't know how to avoid!

If you have done that you should be able to run the command like this:

```
snakemake -s /home3/redwards/GitHubs/EdwardsLab/assembly/mummerplot.snakefile > snakemake.out 2>&1
```


