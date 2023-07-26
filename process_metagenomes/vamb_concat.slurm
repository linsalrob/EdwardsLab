#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH -o vamb_concat-%j.out
#SBATCH -e vamb_concat-%j.err

mkdir -p vamb
python /home/edwa0468/GitHubs/EdwardsLab/process_EK_metagenomes/vamb_concatenate.py vamb/contigs.fna.gz megahit/\*/output/final.contigs.fa

