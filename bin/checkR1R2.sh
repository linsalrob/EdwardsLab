#!/bin/bash


if [[ $# != 1 ]]; then echo "$0  <directory of fastq files>" >&2; exit; fi

# Check that we have an R2 for every R1 and vice-versa

for R1 in $(find $1 -name \*R1\*); do
	R2=${R1/R1/R2};
	if [[ ! -e $R2 ]]; then echo "Not found: $R2 for associated R1: $R1" >&2; fi
done


for R2 in $(find $1 -name \*R2\*); do
	R1=${R2/R2/R1};
	if [[ ! -e $R1 ]]; then echo "Not found: $R1 for associated R2: $R2" >&2; fi
done

