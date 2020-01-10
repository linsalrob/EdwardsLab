use strict;

my $f1 = shift || die "file 1?";
my $f2 = shift || die "file 2?";

my %data; my %first;
open(IN, $f1) || die "cant open $f1";
while (<IN>) {
	$data{$_}=1;
	my @a=split /\t/;
	$first{$a[0]}=$_;
}
close IN;

open(IN, $f2) || die "cant open $f2";
while (<IN>) {
	if (!$data{$_}) {
		my @a=split /\t/;
		if ($first{$a[0]}) {
			print "\n$first{$a[0]}$_\n";
		}
		else {
			print STDERR "NONE: $_";
		}
	}
}
