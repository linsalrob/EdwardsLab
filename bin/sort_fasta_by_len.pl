#!/usr/bin/env perl

=pod

Just sort a fasta file by the length of the sequences.

Default is longest -> shortest.

=cut

use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;
my $rob = new Rob;

my %opts;
getopts('f:vr', \%opts);
unless ($opts{f}) {
	die <<EOF;
	$0
	-f fasta file to parse (required)
	-r sort shortest to longest (default is longest first)
	-v verbose output
EOF
}

my $fa = $rob->read_fasta($opts{f});
my @keys;
if ($opts{r}) {@keys = sort {length($fa->{$a}) <=> length($fa->{$b})} keys %$fa}
else {@keys = sort {length($fa->{$b}) <=> length($fa->{$a})} keys %$fa}

map {print ">$_\n", $fa->{$_}, "\n"} @keys;




