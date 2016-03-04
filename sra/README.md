# Some analysis of kmers, percent prok, and percent phage in the SRA data.

Kate has written PARTIE tha identifies these features in the SRA data, and this is some python code to visualize the data.

The two data files are [run_accession-experiment_lib.tsv.gz](run_accession-experiment_lib.tsv.gz) which has the SRA run accession number and the library type which should be one of:
* AMPLICON
* Bisulfite-Seq
* ChIP-Seq
* CLONE
* CLONEEND
* CTS
* EST
* FINISHING
* FL-cDNA
* MBD-Seq
* miRNA-Seq
* OTHER
* POOLCLONE
* RAD-Seq
* RNA-Seq
* Synthetic-Long-Read
* Tn-Seq
* WCS
* WGA
* WGS
* WXS

And [SRA.partie.tsv.gz](SRA.partie.tsv.gz) which is Kate's original data output.

* filter.py just reads the SRA runs and filters them based on the kmer count or the 16S count that you specify. (The default values are not to filter anything).
* plot_partie_3d.py makes a 3D plot of the percent prokaryote, percent 16S genes, and kmer frequncy
* plot_partie_boxes.py makes a box/whisker plot of the data.