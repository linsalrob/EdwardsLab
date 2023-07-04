

mkdir --parents megahit
for R1 in $(cat R1_reads.txt); do
	R2=${R1/R1/R2};
	FILEEND="_R1_001.fastq.gz";
	O=${R1/$FILEEND/};


	if [[ ! -e megahit/$O ]]; then
		sbatch ~/GitHubs/EdwardsLab/process_EK_metagenomes/megahit.slurm no_sharks/$R1 no_sharks/$R2 megahit/$O
	fi;
done
