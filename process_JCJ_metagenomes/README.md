# Process the CF metagenomes


0. Download the metagenomes

We renamed them to Jess' names. This might not be the best idea, since it makes it harder to backtack when discussing with SAGC, but it makes it easier for subsequent discussions.

Several of the steps will use the file `R1_reads.txt` to know which reads we have. Of course, `snakemake` doesn't need this, but at the moment I'm being very conservative and running each command separately via `slurm`. In addition, this allows me to use `$BGFS` which I don't know how to use with `snakemake`!

```
find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt
NUM_R1_READS=$(wc -l R1_reads.txt | cut -f 1 -d ' ')
echo "There are $NUM_R1_READS R1 reads"
if [[ $(find fastq -name \*R2\* | awk 's++{} END {print s}') != $NUM_R1_READS ]]; then echo "There are a different number of R1 and R2 reads"; fi
```


**Does the length of this match the number of samples we expect?**


1. Remove adapters, exclude reads with > 1 `N`. and filter to reads longer than 100 bp

We use `fastp` for this step, and the adapters are from Illumina. SAGC should be doing this, but they also use `cutadapt` so this gives us belt and braces approaches.


```
sbatch --array=1-$NUM_R1_READS:1 ~/GitHubs/EdwardsLab/process_JCJ_metagenomes/fastp.slurm
```


2. Filter the `fastp` trimmed sequences against the human genome


Next, we want to remove any human genomes from the data. We map the reads to the human genome, and then use our [samtools flags](https://edwards.flinders.edu.au/command-line-deconseq/) to filter human and non-human sequences.

Note that we are using the NCBI human genomes for pipelines for this 

```
sbatch --array=1-$NUM_R1_READS:1 ~/GitHubs/EdwardsLab/process_JCJ_metagenomes/humans.slurm
```

3. Convert the `fastq` sequences to `fasta` sequences. 

Note that `mmseqs easy-taxonomy` requires a fasta file, so here we convert the fastq to fasta files. We also take advantage of this to either add /1 /2 to the R1/R2 reads respectively. It helps with subsequence parsing!!

```
sbatch ~/GitHubs/EdwardsLab/process_JCJ_metagenomes/fastq2fasta.slurm no_human
```

7. Run `mmseqs` against the UniRef50 database

Note that this is  bash script that submits slurm jobs

```
bash ~/GitHubs/EdwardsLab/process_JCJ_metagenomes/mmseqs_easy_taxonomy_submit.sh UniRef50 mmseqs fasta

```









# Assembly and Binning

9. Assembly

I'm using `megahit` for this, just because it works and finishes in a reasonable time. An easy megahit submit script takes care of the hard work:

```
bash /home/edwa0468/GitHubs/EdwardsLab/process_EK_metagenomes/megahit_submit.sh
```

This creates separate directories in the `megahit` directory with each assembly.


