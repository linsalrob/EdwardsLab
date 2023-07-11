#!/bin/bash
#SBATCH --job-name=CountMMSeqs
#SBATCH --time=0-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -o count_mmseqs-%j.out
#SBATCH -e count_mmseqs-%j.err
#SBATCH --partition=short

eval "$(conda shell.bash hook)"
conda activate bioinformatics

perl ~/GitHubs/EdwardsLab/process_EK_metagenomes/count_mmseqs.pl -d mmseqs 

