# Growth in Australian Metagenomes Over Time

How many Australian metagenomes are there?

How frequently are new metagenomes being created or released?

We're going to take a look at the [SRA](https://www.ncbi.nlm.nih.gov/sra) to answer this question.

# Get the data!

We use Google's bigquery to get the data from the SRA table.

Our query is:

```
select * from `nih-sra-datastore.sra.metadata`  where geo_loc_name_country_calc  = 'Australia' 
and (librarysource = "METAGENOMIC" or librarysource = 'METATRANSCRIPTOMIC' or organism like "%microbiom%" OR organism like "%metagenom%");
```

This gets all the results and you can download them as a `.json` file. 

Here is the query:


```
select releasedate,mbytes from `nih-sra-datastore.sra.metadata`  where geo_loc_name_country_calc  = 'Australia' and (librarysource = "METAGENOMIC" or librarysource = 'METATRANSCRIPTOMIC' or organism like "%microbiom%" OR organism like "%metagenom%");
```

You can click on `SAVE RESULTS` to save this as a CSV file to Google Drive, and I suggest renaming the file to something you can find again

# Plot the data

Go to [Google Colab](https://colab.research.google.com/) and create a new notebook. Import the usual panda/seaborn stuff:

```


