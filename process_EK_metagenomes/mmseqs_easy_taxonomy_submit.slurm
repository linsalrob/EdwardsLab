#!/bin/bash
#SBATCH --job-name=mmseqs_submit
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -o mmseqs_submit-%j.out
#SBATCH -e mmseqs_submit-%j.err

# this is so we can submit it as a slurm job :)

bash process_EK_metagenomes/mmseqs_easy_taxonomy_submit.sh UniRef50 mmseqs fasta

