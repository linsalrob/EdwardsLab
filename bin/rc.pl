#!/usr/bin/perl -w
#

use strict;
use Rob;
my $r = new Rob;

my $f = shift || die "Fasta file to reverse complement?";
my $fa = $r->read_fasta($f);
foreach my $id (keys %$fa) {
	print ">$id\n", $r->rc($fa->{$id}), "\n";
}
