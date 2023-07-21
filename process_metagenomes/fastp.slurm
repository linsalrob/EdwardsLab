#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o fastp-%A_%a.out
#SBATCH -e fastp-%A_%a.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics


if [[ ! -e R1_reads.txt ]]; then 
	echo "Please make a file with the R1 reads using this command:" >&2
	echo "find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt" >&2;
	exit 2;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

mkdir -p $BGFS/output
mkdir -p $BGFS/fastq
cp /home/edwa0468/GitHubs/fast-adapter-trimming/adapters/IlluminaAdapters.fa $BGFS

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/R1/R2}
cp $SOURCE/$R1 $SOURCE/$R2 $BGFS/fastq

echo "fastp -n 1 -l 100 -i $BGFS/fastq/$R1 -I $BGFS/fastq/$R2 -o $BGFS/output/$R1 -O $BGFS/output/$R2 --adapter_fasta $BGFS/IlluminaAdapters.fa " >&2
fastp -n 1 -l 100 -i $BGFS/fastq/$R1 -I $BGFS/fastq/$R2 -o $BGFS/output/$R1 -O $BGFS/output/$R2 --adapter_fasta $BGFS/IlluminaAdapters.fa 

mkdir -p fastq_fastp
mv $BGFS/output/$R1 $BGFS/output/$R2  fastq_fastp
