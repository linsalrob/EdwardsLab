#!/bin/bash
#SBATCH --job-name=count_sharks
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o count_sharks-%j.out
#SBATCH -e count_sharks-%j.err

for F in no_sharks/*; do echo -ne "$F\t"; gunzip -c $F  | perl -ne 'if (/^@/) {$s{$_}++} END {print scalar(keys %s), "\n"}'; done > count_nosharks.txt &
for F in    sharks/*; do echo -ne "$F\t"; gunzip -c $F  | perl -ne 'if (/^@/) {$s{$_}++} END {print scalar(keys %s), "\n"}'; done > count_sharks.txt &

wait;

