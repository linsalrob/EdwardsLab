#!/bin/bash
#SBATCH --job-name=Human
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o human-%A_%a.out
#SBATCH -e human-%A_%a.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics


if [[ ! -e R1_reads.txt ]]; then
	echo "Please make a file with the R1 reads using this command:" >&2
	echo "find fastq -name \*R1\* > R1_reads.txt" >&2;
	exit 2;
fi

SOURCE=fastq_fastp
HUMAN=/home/edwa0468/Databases/human/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

mkdir -p $BGFS/output $BGFS/human $BGFS/no_human
mkdir -p $BGFS/fastq
cp $HUMAN $BGFS/human.fna.gz

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/R1/R2}

echo "Processing $R1 and $R2" >&2;
echo "Processing $R1 and $R2";

cp $SOURCE/$R1 $SOURCE/$R2 $BGFS/fastq

# 1651490_20171010_S_R1.fastq.gz
FILEEND="_R1.fastq.gz"
O=${R1/$FILEEND/.bam}
if [[ $O == $R1 ]]; then
	echo "ERROR: OUTPUT $O and INPUT $R1 are the same. Are you sure the R1 files end $FILEEND?" >&2;
	exit 1;
fi

minimap2 -t 16 --split-prefix=$BGFS/tmp$$ -a -xsr $BGFS/human.fna.gz $BGFS/fastq/$R1 $BGFS/fastq/$R2 | samtools view -bh | samtools sort -o $BGFS/output/$O -
samtools index $BGFS/output/$O



echo "From $O"
echo "R1 matching human genome:"
samtools view -F 3588 -f 65 $BGFS/output/$O | cut -f 1 | awk '!s[$1]++' | wc -l
samtools fastq -F 3588 -f 65 $BGFS/output/$O | gzip -c > $BGFS/human/$R1
echo "R2 matching human genome:"
samtools view -F 3588 -f 129 $BGFS/output/$O | cut -f 1 | awk '!s[$1]++' | wc -l
samtools fastq -F 3588 -f 129 $BGFS/output/$O | gzip -c > $BGFS/human/$R2




samtools fastq -F 3584 -f 77 $BGFS/output/$O | gzip -c > $BGFS/no_human/$R1
samtools fastq -F 3584 -f 141 $BGFS/output/$O | gzip -c > $BGFS/no_human/$R2

# We don't currently save singleton reads, but these are they:
# samtools fastq -f 4 -F 1 $SAMPLE.bam > $SAMPLE.unmapped.single.fastq

mkdir -p human no_human
cp -r  $BGFS/human/* human/
cp -r  $BGFS/no_human/* no_human/
