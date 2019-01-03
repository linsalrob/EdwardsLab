# crAssphage

There are several different things in here, but one key is the generator of the files to submit to the NCBI.

# NCBI Submission

We start with a single bioproject ID: PRJNA510571

Then we create biosamples and and id files:

```
python3 ~/GitHubs/EdwardsLab/crAssphage/NCBI_submission.py -d csv -n biosamples.tsv -i id_samples.tsv
```

Then, we submit those BioSamples to NCBI using their website. This gives us two files back via email per submission, one called Objects.txt and one called BioSampleObjects.txt. These files include the biosample names. Since we have three files, I put them in the directory `BioSampleSubmission` and called them starting with the first biosample name: `BioSampleSubmission/SAMN10656826_BioSampleObjects.txt`, `BioSampleSubmission/SAMN10657727_BioSampleObjects.txt`, and `BioSampleSubmission/SAMN10658653_BioSampleObjects.txt`

Now we need to add that data, all the sequences, and everything, and make a single tab separated file that we can use to parse out the sequences:

```
python3 ~/EdwardsLab/crAssphage/NCBI_add_biosample_to_tsv.py -b BioSampleSubmission/SAMN10656826_BioSampleObjects.txt -b BioSampleSubmission/SAMN10657727_BioSampleObjects.txt -b BioSampleSubmissio /SAMN10658653_BioSampleObjects.txt -i BioSampleSubmission/id_samples.tsv -c csv -e all_data.tsv  > all_data_biosamples.tsv
```

Finally, we convert those to SRA metadata files that we can upload to the SRA site, and fasta files of the DNA sequences, one file per biosample.

```
python3 ~/GitHubs/EdwardsLab/crAssphage/NCBI_SRA_Submission.py -f all_data_biosamples.tsv -o sra_metadata.txt -d SRA_Submissions
```

This makes a lot of fasta files in 3 directories, called 0, 1, and 2. The easiest way to upload these was to put the fasta files in an S3 bucket without the subdirectories, and have NCBI download them from there.

Good luck!
