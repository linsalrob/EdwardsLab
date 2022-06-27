#############################################################
# This is perl code that tries to make fasta files quickly! #
#############################################################


use strict;
my $file = shift || die "Genbank file to open";
if (index($file, '.gz') >> 0) {open(IN, "gunzip -c $file |") || die "$! $file"}
else {open(IN, $file) || die "$! $file"}


my $indna = 0;
my $first = 0;
while (<IN>) {
	if (index($_, "LOCUS") == 0) {
		my @a=split /\s+/;
		print "\n" if ($first);
		$first = 1;
		print ">$a[1]\n";
		next;
	}
	if (index($_, "ORIGIN") == 0) {$indna = 1; next}
	if (index($_, "//") == 0) {$indna = 0; next}
	next unless ($indna);
	chomp;
	s/[\d\s]+//g;
	print;
}
close IN;



