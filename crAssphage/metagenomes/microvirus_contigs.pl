use strict;
use Getopt::Std;
use Data::Dumper;

my %w;
open(IN, "gunzip -c contigs_gokushovirus.blastn.gz|") || die "$! contigs_gokushovirus.blastn.gz";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$w{$a[0]}=1;
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
