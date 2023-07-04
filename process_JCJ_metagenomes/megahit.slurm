#!/bin/bash
#SBATCH --job-name=MegaHit
#SBATCH --time=7-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH -o megahit-%j.out
#SBATCH -e megahit-%j.err

# NOTE: 
# - You need the --qos option specified
# - Maxiumum runtime: 7 days
# - 

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

if [[ $# != 3 ]]; then
	echo "$0 <R1 reads> <R2 reads> <output directory>" >&2;
	exit 2;
fi

eval "$(conda shell.bash hook)"
conda activate bioinformatics

mkdir $BGFS/input 
cp -t $BGFS/input $1 $2
R1=$BGFS/input/`basename $1`;
R2=$BGFS/input/`basename $2`;
O=$BGFS/output

megahit -1 $R1 -2 $R2 -o $O

megahit_toolkit contig2fastg 119 $BGFS/output/final.contigs.fa > $BGFS/output/final.graph.fastg; 
fastg2gfa $BGFS/output/final.graph.fastg > $BGFS/output/final.gfa

mkdir -p `dirname $3`

rsync -av $O $3
