#!/bin/bash
#SBATCH --job-name=CmprssOutput
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -o cmprss-%j.out
#SBATCH -e cmprss-%j.err

set -euo pipefail
# here the brackets are required otherwise -exec will only do the last one
find ReadAnnotations/ \( -name "*.tsv" -o -name "*.m8" -o -name "*.xls" \) -exec pigz {} \;

