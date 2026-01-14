#!/bin/bash
#SBATCH --job-name=aa23di
#SBATCH --time=0-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --mem=20G
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

eval "$(conda shell.bash hook)"
conda activate orf_entropy

MODEL="Rostlab/ProstT5"

if [[ $# -gt 0 ]]; then
    MODEL="$1"
fi

python aa23di.py --model "$MODEL" 





