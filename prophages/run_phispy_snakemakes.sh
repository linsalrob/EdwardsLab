#  shell script so I can run all the snakemakes!

WD=$PWD
cd phispy_metrics
echo "Running phispy in phispy_metrics"
snakemake -s phispy_metrics.snakefile -j 12
snakemake -s phispy_no_metrics.snakefile -j 12
python3 summarize.py
cd $WD

cd phispy_tests
echo "Running phispy in phispy_tests"
snakemake -s /home3/redwards/GitHubs/EdwardsLab/prophages/phispy_training_vs_test.snakefile -j 12
cd $WD


cd phispy_training_set
echo "Running phispy in phispy_training_set";
snakemake -s /home3/redwards/GitHubs/EdwardsLab/prophages/phispy_with_training.snakefile -j 12
cd $WD

cd phispy_phage_genes
echo "Running phispy in phispy_phage_genes"
snakemake -s /home3/redwards/GitHubs/EdwardsLab/prophages/phispy_phage_genes.snakefile -j 12
cd $WD


cd PhiSpy_SN
echo "Running phispy in phispy_SN"
snakemake -s /home3/redwards/GitHubs/EdwardsLab/snakemake/phispy.snakefile -j 12
cd $WD

