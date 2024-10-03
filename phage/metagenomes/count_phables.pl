use strict;
use Getopt::Std;
use Data::Dumper;


open(IN, "mv_sequences.txt") || die "$! mv_sequences.txt";
my %class;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$class{$a[0]}=$a[1];
}
close IN;

my $count;
my $header;
open(IN, "sample_genome_read_counts.tsv") || die "$! sample_genome_read_counts.tsv";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	unless ($header) {$header=\@a; next}
	if ($class{$a[0]}) {
		map {$count->{$class{$a[0]}}->{$header->[$_]}=1 if ($a[$_])} (1..$#a);
	}
}
close IN;


foreach my $c (keys %$count) {
	foreach my $s (keys %{$count->{$c}}) {
		print "$c\t$s\n";
	}
}

