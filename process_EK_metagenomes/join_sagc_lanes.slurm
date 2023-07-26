#!/bin/bash
#SBATCH --job-name=JoinSAGCLanes
#SBATCH --time=0-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -o joinlanes-%j.out
#SBATCH -e joinlanes-%j.err
#SBATCH --partition=short

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics


mkdir --parents COMBINED $BGFS/COMBINED/

for L1R1 in $(find L01 -name \*R1\* -printf "%f\n"); do
	L1R2=${L1R1/R1/R2}
	
	L2R1=${L1R1/L01/L02}
	L2R2=${L1R2/L01/L02}

	L3R1=${L1R1/L01/L03}
	L3R2=${L1R2/L01/L03}

	L4R1=${L1R1/L01/L04}
	L4R2=${L1R2/L01/L04}

	if [[ ! -e L01/$L1R1 ]]; then echo "L01/$L1R1 not found" >&2; exit; fi
	if [[ ! -e L01/$L1R2 ]]; then echo "L01/$L1R2 not found" >&2; exit; fi

	if [[ ! -e L02/$L2R1 ]]; then echo "L02/$L2R1 not found" >&2; exit; fi
	if [[ ! -e L02/$L2R2 ]]; then echo "L02/$L2R2 not found" >&2; exit; fi

	if [[ ! -e L03/$L3R1 ]]; then echo "L03/$L3R1 not found" >&2; exit; fi
	if [[ ! -e L03/$L3R2 ]]; then echo "L03/$L3R2 not found" >&2; exit; fi

	if [[ ! -e L04/$L4R1 ]]; then echo "L04/$L4R1 not found" >&2; exit; fi
	if [[ ! -e L04/$L4R2 ]]; then echo "L04/$L4R2 not found" >&2; exit; fi

	OR1=${L1R1/_L01/}
	OR2=${L1R2/_L01/}

	echo "cat L01/$L1R1 L02/$L2R1 L03/$L3R1 L04/$L4R1 > $BGFS/COMBINED/$OR1" >&2;
	echo "cat L01/$L1R2 L02/$L2R2 L03/$L3R2 L04/$L4R2 > $BGFS/COMBINED/$OR2" >&2;

	cat L01/$L1R1 L02/$L2R1 L03/$L3R1 L04/$L4R1 > $BGFS/COMBINED/$OR1
	cat L01/$L1R2 L02/$L2R2 L03/$L3R2 L04/$L4R2 > $BGFS/COMBINED/$OR2
	
	mv $BGFS/COMBINED/* COMBINED/
done

