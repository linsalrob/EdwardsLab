use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;



open(IN, "mv_samples.txt") || die "$! mv_samples.txt";
my %sample;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$sample{$a[0]}=$a[1];
}
close IN;

open(IN, "mv_sequences.txt") || die "$! mv_sequences.txt";
my %contigs;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	if ($sample{$a[0]}) {print "$_\t$sample{$a[0]}\n"}
	elsif (/^contig/) {$contigs{$a[0]}=$_}
	else {print STDERR "Huh? $_\n"}
}
close IN;

open(IN, "gunzip -c contig_count_table.tsv.gz|") || die "$! contig_count_table.tsv.gz";
my $contsamps;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	if ($a[3] > 0) {$contsamps->{$a[1]}->{$a[0]}=1}
}
close IN;

foreach my $c (keys %contigs) {
	print "$contigs{$c}\t";
	if ($contsamps->{$c}) {
		print scalar(keys %{$contsamps->{$c}}), "\n";
	} else {
		print "UNKNOWN\n";
	}
}


