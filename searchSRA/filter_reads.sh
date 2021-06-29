#!/bin/bash

#####################################################################################
#                                                                                   #
# Filter reads, written by Rob Edwards, 28/6/21                                     #
#                                                                                   #
# Note: This is also available as a snakemake pipeline, but that really craps out   #
#       since we have a very large number of small files to process. Find is a lot  #
#       faster!                                                                     #
#                                                                                   #
# To run this, use                                                                  #
#                                                                                   #
#                                                                                   #
#                                                                                   #
#                                                                                   #
#####################################################################################


set -euo pipefail

results="results"
zipfile="results.zip"
mapq=3
length=50
verbose=n
version=0.1
outdir='./'


usage=$(cat <<-EOF
$0
Version $version

Please provide one of either:
-z --zip      The path to the file (usually called results.zip) that you downloaded from SearchSRA
-r --results  The directory with the uncompressed results (e.g. if you have extracted results.zip)
-o --outdir   The directory to write the results (default: $outdir)

-l --length   Minimum alignment length for the sequences to keep the alignment. Default: $length
-m --mapq     Minimum MAPQ (mapping quality) score to keep the alignment. Default $mapq

-v --verbose  More output

-h --help     Print this message and exit

EOF
)

# make sure we have the accessory script!
WD=$(dirname $0)
if [[ ! -e $WD/merge_counts_abstracts.py ]]; then echo "FATAL: $WD/merge_counts_abstracts.py was not found. Can not merge the data"; exit 2; fi
if [[ ! -e $WD/searchSRA_abstracts.tsv.gz ]]; then echo "FATAL: $WD/searchSRA_abstracts.tsv.gz was not found. Can not merge with abstracts"; exit 2; fi


echo -e "Welcome to filter_reads.sh version $version.\nPart of the SearchSRA toolkit.\nCopyright Rob Edwards, 2021"
echo "Started at "$(date)

OPTIONS=r:z:l:m:vho:
LONGOPTS=results:,zip:,length:,mapq:,verbose,help,outdir

PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
	exit 2
fi

eval set -- "$PARSED"

# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
	-h|--help)
	    echo "$usage";
	    exit 0;
	    ;;
        -z|--zip)
            zipfile="$2"
            shift 2
            ;;
        -m|--mapq)
            mapq="$2"
            shift 2
            ;;
        -l|--length)
            length="$2"
            shift 2
            ;;
        -r|--results)
	    results="$2"
            shift 2
            ;;
        -o|--outdir)
	    outdir="$2"
	    mkdir -p $outdir
            shift 2
            ;;
        -v|--verbose)
            verbose=y
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Crap. Some other error"
            exit 3
            ;;
    esac
done

if [[ $# -ne 0 ]]; then echo "$0, don't know what $1 is"; exit 2; fi


echo -n "Results directory: '$results' "
if [[ -d $results ]]; then
	echo "found"; 
else
 	echo "NOT FOUND";
	echo -n "Zipfile: '$zipfile' "
	if [[ -e $zipfile ]]; then
		echo "found. Uncompressing $zipfile to $results";
		unzip -d $results $zipfile;
 	else 
		echo -e "NOT FOUND\n";
		echo "Sorry, please provide either the results directory if you have already uncompressed it"
		echo "or the path to the results zip file you downloaded from search SRA"
		echo "$usage"
		exit 2;
	fi
fi

echo  "Filtering $results using length $length mapq: $mapq"

# Step 1. Filter the reads
# Step 2. Calculate IDX stats

FQDIR=$outdir/FILTERED_MAPQ${mapq}_LEN${length}
IDXDIR=$outdir/IDXSTATS_MAPQ${mapq}_LEN${length}
filterbam=y
idxbam=y

if [[ -d $FQDIR ]]; then echo "$FQDIR already exists, so not filtering the bam files"; filterbam=n; fi
if [[ -d $IDXDIR ]]; then echo "$IDXDIR already exists, so not indexing the bam files"; idxbam=n; fi

if [[ $filterbam == y ]] || [[ $idxbam == y ]]; then
	if [[ $verbose == y ]]; then echo "Parsing and filtering the bam files. Each . is 1000 files"; fi
	mkdir -p $FQDIR $IDXDIR
	C=0
	for BAMF in $(find $results -name \*bam); do
		C=$((C+1))
		if [[ $verbose == y ]]; then
			if [[ $((C % 10000)) == 0 ]]; then echo -n " ";
			elif [[ $((C % 1000)) == 0 ]]; then echo -n "."; fi
		fi
		BAMOUT=$(echo $BAMF | sed -s 's#^.*/##');
		SAMP=$(echo $BAMOUT | sed -e 's/.bam//');
		if [[ $filterbam == y ]]; then
			samtools view -hq $mapq $BAMF | awk -v l=$length 'length($10) > l || $1 ~ /^@/' | samtools view -bS -o $FQDIR/$BAMOUT 
			samtools index $FQDIR/$BAMOUT;
		fi
		if [[ $idxbam == y ]]; then
			samtools idxstats $FQDIR/$BAMOUT | awk '$3 != 0' | cut -f 1,3 | sed -e "s|^|$SAMP\t|" > $IDXDIR/$SAMP.tsv
		fi
	done;

	if [[ $verbose == y ]]; then echo; fi
fi

# Step 3. Combine everything with the python script
if [[ verbose == y ]]; then echo "Running the python merge script to create reads_per_sample.Q${mapq}.L${length}.tsv and reads_per_project.Q${mapq}.L${length}.tsv"; fi
python3 $WD/merge_counts_abstracts.py -d $IDXDIR -a $WD/searchSRA_abstracts.tsv.gz -r $outdir/reads_per_sample.Q${mapq}.L${length}.tsv -p $outdir/reads_per_project.Q${mapq}.L${length}.tsv
