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
	python /home/edwa0468/GitHubs/EdwardsLab/process_JCJ_metagenomes/vamb_concatenate.py vamb/contigs.fna.gz megahit/*/output/final.contigs.fa 
fi


FILEEND="_R1.fastq.gz"
SOURCE=human
OUTDIR=vamb

mkdir --parents $BGFS/$OUTDIR $BGFS/$OUTDIR/mapped_reads  $BGFS/fastq
mkdir --parents $OUTDIR/mapped_reads

# this is needed
cp $OUTDIR/contigs.fna.gz $BGFS/$OUTDIR
if [[ -e $OUTDIR/contigs.mmi ]]; then
	cp $OUTDIR/contigs.mmi $BGFS/$OUTDIR
else
	minimap2 -I100G -d $BGFS/$OUTDIR/contigs.mmi $BGFS/$OUTDIR/contigs.fna.gz
fi

for R1 in $(sort -R R1_reads.txt); do
	R2=${R1/R1/R2}
	BAM=${R1/$FILEEND/.bam}
	if [[ -e $OUTDIR/mapped_reads/$BAM ]]; then
		cp $OUTDIR/mapped_reads/$BAM $BGFS/$OUTDIR/mapped_reads/$BAM;
	else
		cp $SOURCE/$R1 $SOURCE/$R2 $BGFS/fastq
		minimap2 -t 16 -N 5 -ax sr $BGFS/$OUTDIR/contigs.mmi --split-prefix $BGFS/mmsplit$$ $BGFS/fastq/$R1 $BGFS/fastq/$R2 | samtools view -F 3584 -b --threads 16  | samtools sort -@ 16 -o $BGFS/$OUTDIR/mapped_reads/$BAM -
	fi
	rsync -a $BGFS/$OUTDIR/ $OUTDIR/
done

rsync -a $BGFS/$OUTDIR/ $OUTDIR/

vamb --outdir $OUTDIR/vamb --fasta $BGFS/$OUTDIR/contigs.fna.gz --bamfiles $BGFS/$OUTDIR/mapped_reads/*bam -o C -p 16 --cuda

python ~/GitHubs/EdwardsLab/process_JCJ_metagenomes/vamb_create_fasta.py $BGFS/$OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters.tsv 20000 $OUTDIR/bins

rsync -a $BGFS/$OUTDIR/ $OUTDIR/
