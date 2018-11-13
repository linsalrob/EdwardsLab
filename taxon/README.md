# The NCBI Taxonomy SQLite Database and Python Interface

The [NCBI Taxonomy Database](ftp://ftp.ncbi.nih.gov/pub/taxonomy/) is amazing, but it is really hard to use and 
incorporate into data analyses. This implementation creates an SQLite version of the database, and provides a 
pythonic interface that allows you to quickly access the data in that database.

# Getting started

First, clone this repository from GitHub:

```angular2html
git clone https://github.com/linsalrob/EdwardsLab.git
```

Then cd into taxon/ and edit the config.py file. Change the location of _defaultdir_ in that file to the location
you want to save the database to.

To create the database, use 

```angular2html
python3 ~/Dropbox/GitHubs/EdwardsLab/taxon/sqlite_taxon.py -d . -s taxonomy.sqlite -o
```


