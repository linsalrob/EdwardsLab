#############################################################
# This is perl code that tries to make fasta files quickly! #
#############################################################


use strict;
my $file = shift || die "Genbank file to open";
if (index($file, '.gz') >> 0) {open(IN, "gunzip -c $file |") || die "$! $file"}
else {open(IN, $file) || die "$! $file"}


my $indna = 0;
my $defn;
my $id;
my $seq;
while (<IN>) {
	if (index($_, "LOCUS") == 0) {
		my @a=split /\s+/;
		$id = $a[1];
		next;
	}
	if (index($_, "DEFINITION") == 0) {
		chomp;
		s/DEFINITION\s+//;
		$defn = $_;
		my $line = <IN>;
		while (index($line, " ") == 0) {
			chomp($line);
			$defn .= $line;
			$line = <IN>;
		}
		$defn =~ s/\s+/ /g;
		next;
	}
	if (index($_, "ORIGIN") == 0) {$indna = 1; next}
	if (index($_, "//") == 0) {
		$indna = 0;
		if ($id && $defn && $seq) {
			print ">$id [$defn]\n$seq\n";
		} elsif ($id && $seq) {
			print STDERR "NO DEFN FOR $id\n";
			print ">$id [NONE]\n$seq\n";
		} elsif ($seq) {
			print STDERR "NO ID for $seq in $file\n";
		}
		($id, $defn, $seq)=(undef, undef, undef);
		next
	}
	next unless ($indna);
	chomp;
	s/[\d\s]+//g;
	$seq .= $_;
}
close IN;



