use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;
my %count;
open(IN, "most_abundant.txt") || die "$! most_abundant.txt";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$count{$a[0]}=$a[1];
}
close IN;

open(IN, "virus_contig_annotations.tsv") || die "C$! virus_contig_annotations.tsv";
while (<IN>) {
	my @a=split /\t/;
	if ($a[0] eq "contigID") {splice @a, 1, 0, "Samples"}
	elsif ($count{$a[0]}) {splice @a, 1, 0, $count{$a[0]}}
	else {splice @a, 1, 0, 0}
	print join("\t", @a);
}
close IN;
