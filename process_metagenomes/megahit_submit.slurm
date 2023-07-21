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

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh



mkdir --parents megahit
for R1 in $(cat R1_reads.txt); do
	R2=${R1/R1/R2};
	O=${R1/$FILEEND/};


	if [[ ! -e megahit/$O ]]; then
		sbatch ~/GitHubs/EdwardsLab/process_metagenomes/megahit.slurm $HOSTREMOVED/$R1 $HOSTREMOVED/$R2 megahit/$O
	fi;
done
