

mkdir --parents megahit
for R1 in $(cat R1_reads.txt); do
	R2=${R1/R1/R2};
	FILEEND="_R1_001.fastq.gz";
	O=${R1/$FILEEND/};


	if [[ ! -e megahit/$O ]]; then
		sbatch ~/GitHubs/EdwardsLab/process_JCJ_metagenomes/megahit.slurm no_human/$R1 no_human/$R2 megahit/$O
	fi;
done
