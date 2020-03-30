"""

This is a non-working snake file.

Here are the steps that I used. However, it doesn't work as snake

1. Download the genome summary from PATRIC
curl -Lo genome_summary.tsv ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_summary

2. Identify complete genomes
perl -F"\t" -lane 'print $F[0] if ($F[4] eq "Complete")' genome_summary.tsv > complete_ids.txt

3. Create download_faa.sh

```
ID=$(head -n $SGE_TASK_ID complete_ids.txt | tail -n 1)
curl -Lo faa/${ID}.PATRIC.faa ftp://ftp.patricbrc.org/genomes/${ID}/${ID}.PATRIC.faa
```

4. Submit this to the cluster and whack patric :)

```
mkdir faa
chmod +x ./download_faa.sh
qsub -cwd -o sge_out -e sge_err -t 1-27173:1 ./download_faa.sh
```

5. Convert protein IDs to md5 and write an id map file. By definition
this step also dereplicates at 100%. We also ignore proteins with
stop codons, and by default proteins < 100 amino acids

a. first make combine_all.sh that has the python command:
```
python3 /home3/redwards/GitHubs/EdwardsLab/proteins/protein_md5.py -d faa -i id.map -o proteins.faa
```

b. now submit this to the cluster with a hold_jid for the download jobs:
```
qsub -cwd -o comb.sge.out -e comb.sge.err -hold_jid 20693 ./combine_all.sh
```

6. Combine these with mmseqs or cd-hit

Currently mmseqs2 crashes on tatabox. cd-hit has taken 3 days and is
still running!


mmseqs2 command:
```
mmseqs cluster --threads 32 --cov-mode 0 --min-seq-id 0.7 proteins.db proteins.70.clusters $(mktemp -d -p .)
```

cd-hit command:
```
TBD
```

"""







rule download_genome_list:
    output:
        "genome_summary.txt"
    shell:
        "curl -Lo {output} ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_summary"

rule complete_genomes:
    input:
        "genome_summary.txt"
    output:
        "complete_ids.txt"
    shell:
        """
            perl -F"\\t" -lane 'print $F[0] if ($F[4] eq "Complete")' {input} > {output}
        """


