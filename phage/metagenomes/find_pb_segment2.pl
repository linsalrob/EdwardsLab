=pod

contig_199070 is only present in these samples: '35536', '35613', '35634', '35658', '38046'

Can we find any more contigs that are only present in those samples, and not present elsewhere?

=cut

use strict;

open(IN, "gunzip -c contig_count_table.tsv.gz|") || die "$! contig_count_table.tsv.gz";
my $h = 1;

my %w = (
	35536 => 1,
	35613 => 1,
	35634 => 1,
	35658 => 1,
	38046 => 1
);

my %want;
my %other;

while (<IN>) {
	if ($h) {$h=0; next}
	my @a=split /\t/;
	next unless ($a[3]);
	if ($w{$a[0]}) {
		$want{$a[1]}++;
	} else {
		$other{$a[1]}++;
	}
}
close IN;

print "Contig\tWAnted\tOther samples\n";
foreach my $c (keys %want) {
	if ($want{$c} == 5) {
		print "$c\t$want{$c}\t$other{$c}\n";
	}
}


