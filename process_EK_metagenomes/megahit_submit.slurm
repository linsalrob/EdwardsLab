#!/bin/bash
#SBATCH --job-name=megahit_submit
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -o megahit_submit-%j.out
#SBATCH -e megahit_submit-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics

/home/edwa0468/slurm/stats.sh


mkdir --parents megahit
for R1 in $(cat R1_reads.txt); do
	R2=${R1/R1/R2};
	#FILEEND="_R1_001.fastq.gz";
	FILEEND="_R1.fastq.gz";
	O=${R1/$FILEEND/};


	if [[ ! -e megahit/$O ]]; then
		sbatch ~/GitHubs/EdwardsLab/process_EK_metagenomes/megahit.slurm no_sharks/$R1 no_sharks/$R2 megahit/$O
	fi;
done
