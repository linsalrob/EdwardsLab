#__perl__
#
# Calculate the average quality score of one or more files
#
# Usage: average_quality_scores.pl [-l] <list of qual files>
#
# -l initiates printing all scores for each file otherwise
# a summary is produced

use strict;
use Rob;
my $rob = new Rob;
my ($total, $n)=(0,0);
my ($min, $max)=(10000, 0);

my $printall = 0;

foreach my $f (@ARGV) {
	if ($f eq "-l") {$printall = 1; next}
	my $qu = $rob->read_fasta($f, 1);
	foreach my $id (keys %$qu) {
		my @qual = split /\s+/, $qu->{$id};
		my $t=0;
		map {$t+=$_} @qual;
		my $av = $t/($#qual+1);
		$printall && print "$id\t$av\n";
		$total+=$t;
		$n+=$#qual+1;
		$av > $max ? $max = $av : 1;
		$av < $min ? $min = $av : 1;
	}
}
print "TOTAL: $total NBASES: $n MINIMUM: $min MAXIMUM: $max AVERAGE: ", $total/$n, "\n";

