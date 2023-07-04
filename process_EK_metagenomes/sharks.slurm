#!/bin/bash
#SBATCH --job-name=Sharks
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o sharks-%A_%a.out
#SBATCH -e sharks-%A_%a.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics


if [[ ! -e R1_reads.txt ]]; then
	echo "Please make a file with the R1 reads using this command:" >&2
	echo "find fastq -name \*R1\* > R1_reads.txt" >&2;
	exit 2;
fi

SOURCE=fastq_fastp
mkdir -p $BGFS/output $BGFS/sharks $BGFS/no_sharks
mkdir -p $BGFS/fastq
cp /home/nala0006/Ref-Shark-Genomes/Ref-SharkAll.fna $BGFS

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/R1/R2}
cp $SOURCE/$R1 $SOURCE/$R2 $BGFS/fastq


O=${R1/_R1_001.fastq.gz/.bam}

minimap2 -t 16 --split-prefix=$BGFS/tmp$$ -a -xsr $BGFS/Ref-SharkAll.fna $BGFS/fastq/$R1 $BGFS/fastq/$R2 | samtools view -bh | samtools sort -o $BGFS/output/$O -
samtools index $BGFS/output/$O



echo "From $O"
echo "R1 matching shark genome:"
# samtools view -F 4 -f 65 $BGFS/output/$O | cut -f 1 | awk '!s[$1]++' | wc -l
samtools fastq -F 3588 -f 65 $BGFS/output/$O | gzip -c > $BGFS/sharks/$R1
echo "R2 matching shark genome:"
# samtools view -F 4 -f 129 $BGFS/output/$O | cut -f 1 | awk '!s[$1]++' | wc -l
samtools fastq -F 3588 -f 129 $BGFS/output/$O | gzip -c > $BGFS/sharks/$R2




samtools fastq -F 3584 -f 77 $BGFS/output/$O | gzip -c > $BGFS/no_sharks/$R1
samtools fastq -F 3584 -f 141 $BGFS/output/$O | gzip -c > $BGFS/no_sharks/$R2

# We don't currently save singleton reads, but these are they:
# samtools fastq -f 4 -F 1 $SAMPLE.bam > $SAMPLE.unmapped.single.fastq

mkdir -p sharks no_sharks
cp -r  $BGFS/sharks/* sharks/
cp -r  $BGFS/no_sharks/* no_sharks/
