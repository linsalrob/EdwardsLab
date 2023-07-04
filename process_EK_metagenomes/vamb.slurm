#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o vamb-%j.out
#SBATCH -e vamb-%j.err
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate vamb

mkdir --parents vamb

if [[ ! -e vamb/contigs.fna.gz ]]; then
	python /home/edwa0468/GitHubs/EdwardsLab/process_EK_metagenomes/vamb_concatenate.py vamb/contigs.fna.gz megahit/*/output/final.contigs.fa 
fi

mkdir --parents $BGFS/vamb $BGFS/vamb/mapped_reads  $BGFS/fastq
mkdir --parents vamb/mapped_reads

FILEEND="_R1_001.fastq.gz";
SOURCE=no_sharks

cp vamb/contigs.fna.gz $BGFS/vamb
minimap2 -I100G -d $BGFS/vamb/contigs.mmi $BGFS/vamb/contigs.fna.gz
for R1 in $(sort -R R1_reads.txt); do
	R2=${R1/R1/R2}
	cp $SOURCE/$R1 $SOURCE/$R2 $BGFS/fastq
	BAM=${R1/$FILEEND/.bam}
	minimap2 -t 16 -N 5 -ax sr $BGFS/vamb/contigs.mmi --split-prefix $BGFS/mmsplit$$ $BGFS/fastq/$R1 $BGFS/fastq/$R2 | samtools view -F 3584 -b --threads 16  | samtools sort -o $BGFS/vamb/mapped_reads/$BAM -
done

rsync -a $BGFS/vamb/ vamb/

vamb --outdir vamb/vamb --fasta $BGFS/vamb/contigs.fna.gz --bamfiles $BGFS/vamb/mapped_reads/*bam -o C -p 16 --cuda

python ~/GitHubs/EdwardsLab/process_EK_metagenomes/vamb_create_fasta.py $BGFS/vamb/contigs.fna.gz vamb/vamb/vae_clusters.tsv 20000 vamb/bins

rsync -a $BGFS/vamb/ vamb/
