# JPLACER Classification

This is experimental code written by Rob Edwards and Mike Doane.

### What we need

Before we begin, we need three things:

- a directory of fastq files
- a classification file that has the fastq file name and then an arbitrary number of classifications. Each 
classification should be tab-separated.
- the jplacer output file from [phylosift](https://github.com/gjospin/PhyloSift)

### What we will produce

We will output several files that you can import into [itol](https://itol.embl.de)

- a newick file with just the tree, with the metagenomes integrated into the leaves
- multibar files

## Step one, integrate the leaves into the tree

In this step, we read the jplacer file and integrate the metagenomes into the tree. We specifically place ambiquous
sequences at each point where they could occur. We realize that this is not correct, _sensu stricto_, but for the 
visualization it provides us what we need. 

Use the command:
```
python3 ~redwards/EdwardsLab/jplacer/parse_rename_write.py -j sharks_fish.jplace -o sharks_fish.nwk -l sharks_fish.leaves
```

to parse the jplace file and create (a) the tree for itol (sharks_fish.nwk), and (b) a list of all the leaves
(sharks_fish.leaves) that we will use in subsequent commands.

## Step two, create our classifications

We need a directory with the fastq files, and then a list of classifications associated with those files.
The list can be of arbitrary length but must be separated with tabs. Missing values should be empty.

You should probably organise this from the highest to lowest classifications (i.e. the left most 
column should have the least number of unique entries). For example, our classification looks like this:

| fastq file | class 1 | class 2 |
| --- | --- | --- | 
| ts.fastq | sharks | thresher shark | 
| ws.fastq | sharks | whale shark |
| bl.fastq | fish | blennie |
| fl.fastq | fish | flounder |

This step also requires access to the [SQLite3 taxnomy database](https://github.com/linsalrob/EdwardsLab/tree/master/taxon)

We generate a file that has all the leaves found in the tree, including the reference sequences, and whether they are 
bacteria, archaea, eukarya, or from your metagenomes. We also append all the classification information you provide.

````
python3 ~redwards/GitHubs/EdwardsLab/jplacer/fastq2ids.py -l sharks_stingray.leaves -c ../fastq_classification.tsv -d ../fastq -o sharks_stingray.leaves.labels -v 2> err
````

## Step three, create color strips for different taxonomic levels

We can create colorstrips for e.g. Bacteria, Archaea, Metagenomes or for Sharks and Fish based on the data in our 
leaves.labels file that we have just created. e.g. to make a file for the Kingdom, use:

```
python3 ~redwards/GitHubs/EdwardsLab/jplacer/create_colorstrip.py -f sharks_stingray.leaves.labels -n 2 -l Kingdom -s 1 -o sharks_stingray.kingdom.colorstrip
```

and to create a similar file for shark or fish, use:

```
python3 ~redwards/GitHubs/EdwardsLab/jplacer/create_colorstrip.py -f sharks_stingray.leaves.labels -n 4 -l "Fish/Shark" -s 2 -o sharks_stingray.fish_shark.colorstrip
```

and to create one for each of the fish/shark species we use:

```
python3 ~redwards/GitHubs/EdwardsLab/jplacer/create_colorstrip.py -f sharks_stingray.leaves.labels -n 5 -l Species -s 3 -o sharks_stingray.species.colorstrip
```

## Step four, count the metagenomes at different levels and create multibars



