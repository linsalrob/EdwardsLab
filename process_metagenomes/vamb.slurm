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

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


OUTDIR=vamb
mkdir --parents $OUTDIR

if [[ ! -e $OUTDIR/contigs.fna.gz ]]; then
	python /home/edwa0468/GitHubs/EdwardsLab/process_EK_metagenomes/vamb_concatenate.py $OUTDIR/contigs.fna.gz megahit/*/output/final.contigs.fa 
fi

mkdir --parents $BGFS/$OUTDIR $BGFS/$OUTDIR/mapped_reads  $BGFS/fastq
mkdir --parents $OUTDIR/mapped_reads

cp $OUTDIR/contigs.fna.gz $BGFS/$OUTDIR
minimap2 -I100G -d $BGFS/$OUTDIR/contigs.mmi $BGFS/$OUTDIR/contigs.fna.gz
for R1 in $(sort -R R1_reads.txt); do
	R2=${R1/R1/R2}
	BAM=${R1/$FILEEND/.bam}
	if [[ -e $OUTDIR/mapped_reads/$BAM ]]; then
		cp $OUTDIR/mapped_reads/$BAM $BGFS/$OUTDIR/mapped_reads/$BAM;
	else
		cp $HOSTREMOVED/$R1 $HOSTREMOVED/$R2 $BGFS/fastq
		minimap2 -t 16 -N 5 -ax sr $BGFS/$OUTDIR/contigs.mmi --split-prefix $BGFS/mmsplit$$ $BGFS/fastq/$R1 $BGFS/fastq/$R2 | samtools view -F 3584 -b --threads 16  | samtools sort -@ 16 -o $BGFS/$OUTDIR/mapped_reads/$BAM -
	fi
	rsync -a $BGFS/$OUTDIR/ $OUTDIR/
done

rsync -a $BGFS/$OUTDIR/ $OUTDIR/

vamb -o C -p 16 --cuda --outdir $OUTDIR/vamb --fasta $BGFS/$OUTDIR/contigs.fna.gz --bamfiles $BGFS/$OUTDIR/mapped_reads/*bam 

python ~/GitHubs/EdwardsLab/process_EK_metagenomes/vamb_create_fasta.py $BGFS/$OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters.tsv 20000 $OUTDIR/bins

rsync -a $BGFS/$OUTDIR/ $OUTDIR/
