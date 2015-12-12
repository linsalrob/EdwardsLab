#!/usr/bin/perl -w
#

use strict;

# convert a genbanktable to fasta using either contig/start/end/strand or gene id

use Getopt::Std;
my %opts;
getopts('lg:i', \%opts);
unless (($opts{i} || $opts{l}) && $opts{g}) {
	die <<EOF;
$0
-g genbank table file (required)
-l use location (conting, start, stop) as ID
-i use gene id as ID

either i or l must be specified (but not both)
EOF
}

if ($opts{l} && $opts{i}) {die "Not both -i and -l"}

open(IN, $opts{g}) || die "$! $opts{g}";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	my $header;
	if ($opts{i}) {
		$header = ">$a[6]\n";
	} else {
		my ($b, $e, $strand) =  ($a[7], $a[8], $a[9]);
		if ($strand < 0) {($b, $e) = ($e, $b)}
		$header = ">$a[0]_${b}_${e}\n";
	}
	unless ($header) {print STDERR "Can't construct a header in $_\n"}
	print $header, $a[10], "\n";
}
