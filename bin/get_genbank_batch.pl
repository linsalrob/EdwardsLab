#!/usr/bin/env perl
use strict;
use LWP::Simple;
use Data::Dumper;

unless ($ARGV[0]) {die "$0 <file of ids> <number to get (default = 100)" }
open(IN, $ARGV[0]) || die "Can't open $ARGV[0]";
my @ids;
while (<IN>) {
	chomp;
	push @ids, $_;
}
close IN;
my $n=100;
if ($ARGV[1]) {$n=$ARGV[1]}

my $time=time-10;
my $url='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&retmode=text&rettype=gb&tool=robsgetter&email=rob@salmonella.org&id=';
open(OUT, ">sequences.gbk") || die "Can't write to sequences.gbk";
my $total=0;
while (@ids) {
	my @pieces = splice(@ids, 0, $n);
	$total+=scalar(@pieces);
	my $urlid=$url . join(",", @pieces);
	while (time-$time < 1.5) {sleep 1}
	print STDERR "Getting upto $total\n";
	$time=time;
	print OUT get($urlid);
}
close OUT;

exit 0;

