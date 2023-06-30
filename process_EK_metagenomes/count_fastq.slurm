#!/bin/bash
#SBATCH --job-name=countfastq
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o count_fastq-%j.out
#SBATCH -e count_fastq-%j.err

DIR=$1

if [[ -z $DIR ]]; then 
	echo "sbatch count_fastq.slurm <directory name>" >&2;
	echo "You can also use a single fastq file, I suppose" >&2;
	exit 1;
fi

# This requires count_fastq from my EdwardsLab github to be installed

count_fastq $DIR
