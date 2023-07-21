#!/bin/bash
#SBATCH --job-name=CountSS
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH -o count_ss-%j.out
#SBATCH -e count_ss-%j.err

perl /home/edwa0468/GitHubs/EdwardsLab/process_EK_metagenomes/count_subsystems.pl -d mmseqs
