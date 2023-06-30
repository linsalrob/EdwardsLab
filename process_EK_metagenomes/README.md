# This is yet another stand alone metagneome processing pipeline, starting from fastq files provided by SAGC

1. After downloading the `zip` file from SAGC and uncompressing it, we have four lanes of data, called `L01`, `L02`, `L03`, and `L04`. We merge those into a single file:

```
sbatch join_sagc_lanes.slurm
```

This will create one fastq file for each R1 and R2. Then, we check we have an R1 for each R2 and vice versa:

*Check for missing R1/R2 files*

We can easily check to see if anything is missing. We are going to do this at multiple steps along the way

```
~/GitHubs/EdwardsLab/bin/checkR1R2.sh COMBINED
```

*Count the sequences*

We count the sequences, and keep a list of those counts so we can see what happens at every step.

```
sbatch ~/slurm/count_fastq.slurm COMBINED
```


*Checking sequence counts*


If you have some previously extracted files, and you want to compare the counts for the re-downloaded data, then count the fastq files using `count_fastq.slurm` and then 

```
IFS=$'\n';  for F in $(cat EagleRay/count_fastq-2158191.out); do echo -ne "$F\t"; B=$(echo $F | awk '{print $1}' | sed -e 's/fastq/COMBINED/; s/S34/.*/'); C=$(grep -P $B SAGCQA0428/count_fastq-2158535.out); echo -e "$B\t$C"; done
```

2. Deploy the reads into the right directories.

Create a new directory for each subsample, and deploy the reads. No code is currently given for this, you are on your own. Sorry.


3. Create a list of all the R1 reads.

Several of the steps will use the file `R1_reads.txt` to know which reads we have. Of course, `snakemake` doesn't need this, but at the moment I'm being very conservative and running each command separately via `slurm`. In addition, this allows me to use `$BGFS` which I don't know how to use with `snakemake`!

```
find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt
NUM_R1_READS=$(wc -l R1_reads.txt | cut -f 1 -d ' ')
echo There are $NUM_R1_READS R1 readsA
if [[ $(find fastq -name \*R2\* | awk 's++{} END {print s}') != $NUM_R1_READS ]]; then echo "There are a different number of R1 and R2 reads"; fi
```

**Does the length of this match the number of samples we expect?**

4. Run `fastp`.

Nowadays, SAGC should be using `cutadapt` to remove the adapters (??) but I rerun `fastp`. This also trims out sequences with an `N` and low quality sequences. The command for a single run is


```
fastp -n 1 -l 100 \
	-i fastq/SAGCFN_22_00734_S72_R1_001.fastq.gz -I fastq/SAGCFN_22_00734_S72_R2_001.fastq.gz \
	-o fastq_fastp/SAGCFN_22_00734_S72_R1_001.fastq.gz -O  fastq_fastp/SAGCFN_22_00734_S72_R2_001.fastq.gz \
	--adapter_fasta /home/edwa0468/GitHubs/fast-adapter-trimming/adapters/IlluminaAdapters.fa
```

Note here I drop any read with a single `N` and require a minimum length of 100 bp. Moreover, `fastp` drops reads where there is no mate, by default, so all the reads end up properly paired.


You can run that for all the files listed in `R1_reads.txt` using this array job:

```
sbatch --array=1-$NUM_R1_READS:1 ~/GitHubs/EdwardsLab/process_EK_metagenomes/fastp.slurm
```


5. Filter the `fastp` trimmed sequences against the shark genome


Next, we want to remove all the shark genomes from the data. We map the reads to the shark genome, and then use our [samtools flags](https://edwards.flinders.edu.au/command-line-deconseq/) to filter shark and non-shark sequences.

```
sbatch --array=1-$NUM_R1_READS:1 ~/GitHubs/EdwardsLab/process_EK_metagenomes/sharks.slurm
```


6. Convert the `fastq` sequences to `fasta` sequences. 

Note that `mmseqs easy-taxonomy` requires a fasta file, so here we convert the fastq to fasta files. We also take advantage of this to either add /1 /2 to the R1/R2 reads respectively. It helps with subsequence parsing!!

```
sbatch ~/GitHubs/EdwardsLab/process_EK_metagenomes/fastq2fasta.slurm no_sharks/
```

7. Run `mmseqs` against the UniRef50 database

Note that this is  bash script that submits slurm jobs





