
use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;

# What are the most abundant contigs in terms of numbers of samples
# they are in
#


my $count;
open(IN, "gunzip -c contig_count_table.tsv.gz|") || die "$! contig_count_table.tsv.gz";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	next unless ($a[3]);
	$count->{$a[1]}->{$a[0]} = 1;
}
close IN;

foreach my $c (keys %$count) {
	print "$c\t", scalar(keys %{$count->{$c}}), "\n";
}
