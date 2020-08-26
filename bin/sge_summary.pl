use strict;

open(IN, "qstat |") || die "$! pipe to qstat";
my $job;
my %type;
while (<IN>) {
	chomp;
	next if (/^job-ID/ || /^\-/);
	my @a=split /\s+/;
	$job->{$a[2]}->{$a[4]}++;
	$type{$a[4]}=1;
}

my @t = sort {lc($a) cmp lc($b)} keys %type;
print join("\t", "Job       ", @t), "\n";

foreach my $j (sort {lc($a) cmp lc($b)} keys %{$job}) {
	print $j;
	foreach my $t (@t) {
		print $job->{$j}->{$t} ? "\t". $job->{$j}->{$t} : "\t0";
	}
	print "\n";
}


