# Code to work with _PhyloSift_ and _pplacer_ output files

This is all experimental and was written for Mike Doane to use!


To generate the output files for itol:

1. Run this command to parse the jplace file and create (a) the tree for itol (sharks_fish.nwk), and (b) a list of all the leaves (sharks_fish.leaves)

```
python3 ~redwards/EdwardsLab/jplacer/parse_rename_write.py -j sharks_fish.jplace -o sharks_fish.nwk -l sharks_fish.leaves
```

2. Then run the command to generate the color strips. This reads the leaves file (sharks_fish.leaves) from the previous command, and the fastq files (use multiple -f, one for each fastq file), and creates a color strip for the inner circle that has the domain (sharks_fish.leaves.annotations.txt), a color strip for the outer circle that has the fastq name (sharks_fish.metagenome.annotations.txt), and a simple tab file that has the domains and fastq files so we don't have to process everything again!


```
python3 ~redwards/EdwardsLab/jplacer/generate_color_strip.py -l sharks_fish.leaves -f fish_cat.fastq -f sharks_cat.fastq -i sharks_fish.leaves.annotations.txt -o sharks_fish.metagenome.annotations.txt -r sharks_fish.reads.domains 
python3 ~redwards/EdwardsLab/jplacer/fastq2color_strip.py -l sharks_fish.reads.domains -d fastq/ -o fish_species2.txt -v > out 2> e
```


3. To parse the tree file and create the distance matrix, we use

```
python3 ~redwards/EdwardsLab/jplacer/tree_to_cophenetic_matrix.py -t sharks_fish.nwk > sharks_fish.distance.tsv
```

4. To create the labelled distance matrix we use:

```
perl join.pl > sharks_fish.distance.labelled.tsv
```

5. Then we use count metagenomes to get the counts at different taxonomic levels:

```
python3 ~redwards/EdwardsLab/jplacer/count_metagenomes.py -t sharks_fish.nwk -d sharks_fish.distance.labelled.tsv -r fish -l r_family -p
```

