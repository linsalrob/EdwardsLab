use strict;

open(IN, "all_data_counts.tsv") || die "$! all_data_counts.tsv";
my %reads;
while (<IN>) {
	next if (/R1/);
	chomp;
	my @a=split /\t/;
	$reads{$a[0]}=$a[1] + $a[2];
}
close IN;

my %genera; my $data;
foreach my $sa (keys %reads) {
	open(IN, "gunzip -c TaxonomyWithSharks/${sa}_report.gz |") || die "$! TaxonomyWithSharks/${sa}_report.gz";
	my $inbact = 0;
	while (<IN>) {
		next if (/no rank/);
		chomp;
		my @a=split /\s+/;
		if ($a[3] eq "superkingdom") {$inbact = 0}
		if ($a[3] eq "superkingdom" && $a[5] eq "Bacteria") {$inbact = 1}
		next unless ($inbact);
		next unless ($a[3] eq "genus");
		my $bn = join(" ", @a[5..$#a]);
		$data->{$sa}->{$bn}=($a[1]/$reads{$sa} * 1e6);
		$genera{$bn}=1;
		# print join("\t", $a[1]/$reads{$sa}, $bn), "\n";
	}
	close IN;
}

my @genera = sort {$a cmp $b} keys %genera;

print join("\t", "", @genera), "\n";
foreach my $sa (keys %reads) {
	print $sa;
	foreach my $g (@genera) {
		$data->{$sa}->{$g} ? print "\t", $data->{$sa}->{$g} : print "\t0";
	} 
	print "\n";
}



