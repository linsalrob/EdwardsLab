# Phage protein genera

This data is generated Dec 2017. We have all phage genomes from GenBank, and have compared those to all proteins in nr:

1\. Download everything from GenBank

```
ncftpget  ftp://ftp.ncbi.nlm.nih.gov/genbank/gbphg*
```

2\. Extract the sequences into separate genbank files
```
python3 ~/EdwardsLab/bin/separate_multigenbank.py -f gbphg4.seq.gz -d phages
python3 ~/EdwardsLab/bin/separate_multigenbank.py -f gbphg3.seq.gz -d phages
python3 ~/EdwardsLab/bin/separate_multigenbank.py -f gbphg2.seq.gz -d phages 
python3 ~/EdwardsLab/bin/separate_multigenbank.py -f gbphg1.seq.gz -d phages 
```

3\. Find complete phage genomes
```
mkdir complete_genomes

cd phages
grep genome * | grep DEFINITION | grep  complete | cut -f 1 -d ':' | xargs -i cp {} ../complete_genomes/
```

4\. Now find those with host information
```
cd ../
mkdir host
cd complete_genomes/
grep \/host= * | cut -f 1 -d ':' | xargs -i cp {} ../host/
```

4b\. Post a quick tweet about how many more there have been since last time:
<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">There are now 5,088 complete <a href="https://twitter.com/hashtag/phage?src=hash&amp;ref_src=twsrc%5Etfw">#phage</a> genomes in GenBank and 3,072 of them have a host tag. Up from 820 phages w/ host when <a href="https://twitter.com/BEDutilh?ref_src=twsrc%5Etfw">@BEDutilh</a> et al tried to computationally predict host in 2016 <a href="https://t.co/tdfTiWqC1D">https://t.co/tdfTiWqC1D</a></p>&mdash; Rob Edwards (@linsalrob) <a href="https://twitter.com/linsalrob/status/944378872031612928?ref_src=twsrc%5Etfw">December 23, 2017</a></blockquote>


5\. Convert those genbank files  to flatfiles
```
C=0; for h in $(ls host); do C=$((C+1)); echo $C  $h; o=$(echo $h | sed -e 's/gbk/txt/'); genbank2flatfile.pl host/$h; mv tbl flatfiles/$o; done
history | tail -n 50 | perl -pe 's/\s+\d+\s+//'
```


6\. Convert those flat files to a single fasta files for blast
```
cut -f 6,10 flatfiles/* | sed -e 's/^/>/; s/\t/\n/' > phage_genes.faa
```

This has 253,260 proteins, and so we blast those against nr on anthill using something like this:

```
split_blast_queries_edwards_blastplus.pl -f phage_genes.faa -n 40 -d blastp.nr -p blastp -N phageblastp -db /home2/db/blast/nr/nr -evalue 1e-5  -outfmt 6 'std qlen slen staxids sscinames' -max_target_seqs 1000 -num_threads 16
```

NOTE: That command screws up the order of the arguments, and so you need to munge the files:

```
cd blastp.nr/
perl -i -npe "s/'6' -evalue 1e-5   'std qlen slen staxids sscinames'/'6 std qlen slen staxids sscinames' -evalue 1e-5/" blphage_ge.*.sh
```

but it turns out that is OK, because it allows you to submit to the cluster using whole nodes:

```
rm -rf sge_output/; mkdir sge_output/
qsub -cwd -pe make 16 -t 1-40:1 -o sge_output/ -e sge_output/ ./blphage_ge.submitall.0.sh
```

7\. convert the blast output to a list of taxonomies:

```
python3 ~/EdwardsLab/phage_protein_blast_genera/blast_tax_to_genera.py -v -f blastp.nr/*.blastp > ids_taxa.txt
```

8\. Get the hosts for those viruses:

```
grep \/host= * | sed -e 's/:\s\+\/host="/\t/; s/"//g; s/.gbk//' > phage_host.tsv
```

This was munged together with a list of all bacteria associated with humans to create phage_host_location.txt


9\. Make a list of sites and number of gnera that are matched

```
python3 ~/EdwardsLab/phage_protein_blast_genera//genera_per_phage_protein.py -d flatfiles/ -i ids_taxa.txt -b > gav.out
```


10\. Make a figure of the data:

```
python3 ~/EdwardsLab/phage_protein_blast_genera/num_prots_vs_taxa.py
```





