#!/bin/bash
#SBATCH --job-name=MegaHit
#SBATCH --time=7-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH -o megahit-%A_%a.out
#SBATCH -e megahit-%A_%a.err


## For an example run on Emma's EagleRay data
# 
#  This suggests that we probably don't need 128G// high capacity,
#  So I downgraded it.
#
# Job ID: 2159634
# Cluster: deepthought
# User/Group: edwa0468/staff
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 32
# CPU Utilized: 2-02:38:56
# CPU Efficiency: 91.74% of 2-07:12:32 core-walltime
# Job Wall-clock time: 01:43:31
# Memory Utilized: 232.00 KB
# Memory Efficiency: 0.00% of 128.00 GB
# 
# Parameters were
# #SBATCH --mem=128G
# #SBATCH --partition=high-capacity
# #SBATCH --qos=hc-concurrent-jobs
#
#

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
		mkdir $BGFS/input 
		cp -t $BGFS/input $HOSTREMOVED/$R1 $HOSTREMOVED/$R2
		TMPO=$BGFS/$O

		megahit -1 $BGFS/input/$R1 -2 $BGFS/input/$R2 -o $TMPO

		megahit_toolkit contig2fastg 119 $TMPO/final.contigs.fa > $TMPO/final.graph.fastg; 
		fastg2gfa $TMPO/final.graph.fastg > $TMPO/final.gfa


		rsync -av $TMPO megahit/$O
	fi
done
