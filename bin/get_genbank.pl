#!/usr/bin/perl
use strict;
use LWP::Simple;
use Data::Dumper;

my $id = shift || die "$0 <genbank id>";


open(OUT, ">${id}_sequences.gbk") || die "Can't write to ${id}_sequences.gbk";
my $url='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&retmode=text&rettype=gb&id=' . $id;

print OUT get($url);
print "$url\n";
exit 0;

