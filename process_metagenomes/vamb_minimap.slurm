#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o map_reads-%A_%a.out
#SBATCH -e map_reads-%A_%a.err

# Just the minimap section of VAMB which does not require GPUs and we can do in parallel

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
	echo "Please create the contigs with vamb_concatenate" >&2;
	exit 2;
	#python /home/edwa0468/GitHubs/EdwardsLab/process_EK_metagenomes/vamb_concatenate.py $OUTDIR/contigs.fna.gz megahit/*/output/final.contigs.fa 
fi


R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/R1/R2}
BAM=${R1/$FILEEND/.bam}

if [[ ! -e $OUTDIR/mapped_reads/$BAM ]]; then
	mkdir --parents $BGFS/$OUTDIR $BGFS/$OUTDIR/mapped_reads  $BGFS/fastq
	mkdir --parents $OUTDIR/mapped_reads

	if [[ -e $OUTDIR/contigs.mmi ]]; then
		cp $OUTDIR/contigs.mmi $BGFS/$OUTDIR/contigs.mmi
	else
		cp $OUTDIR/contigs.fna.gz $BGFS/$OUTDIR
		minimap2 -I100G -d $BGFS/$OUTDIR/contigs.mmi $BGFS/$OUTDIR/contigs.fna.gz
	fi

	cp $HOSTREMOVED/$R1 $HOSTREMOVED/$R2 $BGFS/fastq
	minimap2 -t 16 -N 5 -ax sr $BGFS/$OUTDIR/contigs.mmi --split-prefix $BGFS/mmsplit$$ $BGFS/fastq/$R1 $BGFS/fastq/$R2 | samtools view -F 3584 -b --threads 16  | samtools sort -@ 16 -o $BGFS/$OUTDIR/mapped_reads/$BAM -
	rsync -a $BGFS/$OUTDIR/ $OUTDIR/
fi
