if [ -z $1 ]; then 
	echo "$0 <gfa file>  >  <fasta file>"
	exit $E_BADARGS
fi

awk -v id=$1 '/^S/{print ">"id"_"$2"\n"$3}' $1
