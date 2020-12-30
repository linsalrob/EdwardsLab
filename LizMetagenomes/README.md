# Processing Metagenomes

Rob Edwards, Aug-Dec 2020

## What is this?

This is a [snakemake](https://snakemake.readthedocs.io/) pipeline for processing metagenomes written by Rob. It should be trivial to run, but most likely won't be.

## Getting started

You will need these files:
   - `process_shark_metagenomes.snakefile` is the snakefile that does the work
   - `liz_metagenomes_env.yaml` is the definition file for all the code we need to use

You should also set up a [conda profile](#setting-your-conda-profile) for anthill or rambox or whereever you work.









## Setting your conda profile

You can set up a profile that snakemake will use automatically. I recommend setting two different profiles, one for anthill and one for rambox. Anthill uses [SGE](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) to submit jobs to the cluster, but rambox is just the command line.

To make these profiles, start by making the two directories:

```
mkdir -p ~/.config/snakemake/anthill ~/.config/snakemake/rambox
```

Then we will make a config file in each. Here are the contents:

In **~/.config/snakemake/anthill/config.yaml** put:

```
jobs: 600
use-conda: True
conda-frontend: mamba
default-resources: [cpus=1, mem_mb=2000]
keep-going: True

# sge settings
cluster: "qsub -cwd -o sge_out -e sge_err -pe smp {resources.cpus} -V "
latency-wait: 60
local-cores: 6
```

In **~/.config/snakemake/rambox/config.yaml** put:

```
jobs: 600
use-conda: True
conda-frontend: mamba
default-resources: [cpus=1, mem_mb=2000]
keep-going: True

# sge settings
local-cores: 24
```


Now, if you are running the code on anthill you can use `--profile anthill` and if you are running the snakemake code on rambox you can use `--profile rambox`
