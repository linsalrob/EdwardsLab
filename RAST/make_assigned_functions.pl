#!/usr/bin/perl 
#

use strict;

my $dir = shift || die "seed directory";
my $odir = shift;
unless ($odir) {$odir = $dir}

if (-e "$dir/assigned_functions") {
	print STDERR "Backing up assigned_functions\n";
	`cp -f $dir/assigned_functions $dir/assigned_functions.bak`;
}

my %fn;
open(IN, "$dir/proposed_non_ff_functions") ||die "$! proposed_non_ff_functions";
while(<IN>) {
	chomp;
	my @a=split /\t/;
	$fn{$a[0]}=$a[1];
}
close IN;

open(IN, "$dir/proposed_functions") ||die "$! proposed_functions";
while(<IN>) {
	chomp;
	my @a=split /\t/;
	$fn{$a[0]}=$a[1];
}
close IN;



open(OUT, ">$odir/assigned_functions") || die "$! assigned_functions";
map {print OUT "$_\t$fn{$_}\n"} sort {$a cmp $b} keys %fn;
close OUT;


