use strict;
use Getopt::Std;
use Data::Dumper;

my %w;
open(IN, "gunzip -c picobirnaviridae_rdrp.tblastn.gz |") || die "$! picobirnaviridae_rdrp.tblastn.gz";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$w{$a[1]}=1;
}
close IN;

open(IN, "gunzip -c contig_count_table.tsv.gz|") || die "$! contig_count_table.tsv.gz";
my $h = 1;
while (<IN>) {
	if ($h) {print; $h=0; next}
	my @a=split /\t/;
	print if ($w{$a[1]} && $a[3]);
}
close IN;
