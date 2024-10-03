use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;

# the first version was shit because it counted hits to multiple viruses more than once. 
# now we need to know which samples have e.g. Gokushaviridae
#

# get the contig information first since I need to work on the phables information
#

open(IN, "mv_sequences.txt") || die "$! mv_sequences.txt";
my %class;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$class{$a[0]}=$a[1];
}
close IN;

my $count;
open(IN, "gunzip -c contig_count_table.tsv.gz|") || die "$! contig_count_table.tsv.gz";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	next unless ($a[3]);
	if ($class{$a[1]}) {
		$count->{$class{$a[1]}}->{$a[0]} = 1;
	}
}
close IN;

foreach my $c (keys %$count) {
	foreach my $s (keys %{$count->{$c}}) {
		print "$c\t$s\n";
	}
}
