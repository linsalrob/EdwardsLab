#!/bin/bash
#SBATCH --job-name=HostRemoval
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o host-%A_%a.out
#SBATCH -e host-%A_%a.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics


if [[ ! -e R1_reads.txt ]]; then
	echo "Please make a file with the R1 reads using this command:" >&2
	echo "find fastq -name \*R1\* > R1_reads.txt" >&2;
	exit 2;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED, HOSTFILE, etc" >&2
	exit 2;
fi

source DEFINITIONS.sh


QC=fastq_fastp

mkdir -p $BGFS/output $BGFS/$HOST $BGFS/$HOSTREMOVED
mkdir -p $BGFS/fastq
cp $HOSTFILE $BGFS/host.fna.gz

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/R1/R2}

echo "Processing $R1 and $R2" >&2;
echo "Processing $R1 and $R2";

cp $QC/$R1 $QC/$R2 $BGFS/fastq

O=${R1/$FILEEND/.bam}
if [[ $O == $R1 ]]; then
	echo "ERROR: OUTPUT $O and INPUT $R1 are the same. Are you sure the R1 files end $FILEEND?" >&2;
	exit 1;
fi

minimap2 -t 16 --split-prefix=$BGFS/tmp$$ -a -xsr $BGFS/host.fna.gz $BGFS/fastq/$R1 $BGFS/fastq/$R2 | samtools view -bh | samtools sort -o $BGFS/output/$O -
samtools index $BGFS/output/$O


# Reads that match the host
samtools fastq -F 3588 -f 65 $BGFS/output/$O | gzip -c > $BGFS/$HOST/$R1
samtools fastq -F 3588 -f 129 $BGFS/output/$O | gzip -c > $BGFS/$HOST/$R2


# Reads that do not match the host

samtools fastq -F 3584 -f 77 $BGFS/output/$O | gzip -c > $BGFS/$HOSTREMOVED/$R1
samtools fastq -F 3584 -f 141 $BGFS/output/$O | gzip -c > $BGFS/$HOSTREMOVED/$R2

# We don't currently save singleton reads, but these are they:
# samtools fastq -f 4 -F 1 $SAMPLE.bam > $SAMPLE.unmapped.single.fastq

mkdir -p $HOST $HOSTREMOVED
rsync -a $BGFS/$HOST/ $HOST/
rsync -a $BGFS/$HOSTREMOVED/ $HOSTREMOVED/
