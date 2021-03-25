#!/usr/bin/env perl
use strict;
use LWP::UserAgent;
use Data::Dumper;

unless ($ARGV[0]) {die "$0 <file of ids> <number to get (default = 100)" }
open(IN, $ARGV[0]) || die "Can't open $ARGV[0]";
my @ids;
while (<IN>) {
	chomp;
	s/\r//g;
	push @ids, $_;
}
close IN;
my $n=100;
if ($ARGV[1]) {$n=$ARGV[1]}


my $ua = LWP::UserAgent->new;
$ua->env_proxy;

my $time=time-10;
my $url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&retmode=text&rettype=gb&tool=robsgetter&email=rob@salmonella.org&id=';
open(OUT, ">sequences.gbk") || die "Can't write to sequences.gbk";
my $total=0;
while (@ids) {
	my @pieces = splice(@ids, 0, $n);
	$total+=scalar(@pieces);
	my $urlid=$url . join(",", @pieces);
	while (time-$time < 1.5) {sleep 1}
	print STDERR "Getting upto $total\n";
	$time=time;

	my $response = $ua->get($urlid);
	if ($response->is_success) {
		print OUT $response->decoded_content;
	} else {
		print STDERR "ERROR getting upto $total. Response code: ", $response->status_line, "\n";
		print STDERR "$urlid\n";
		exit(-1);
	}


}
close OUT;

exit 0;

