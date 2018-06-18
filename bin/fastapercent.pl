#!/usr/bin/perl -w

# look through a fasta file and figure out what percent is done

use strict;

my ($file, $tag)=@ARGV;
unless ($file && $tag) {die "$0 <fasta file> <fasta tag (without '>')>"}
my ($c, $s)=(0,0);
open(IN, $file) || die "can't open $file";
while (<IN>) {
 next unless (/^>/);
 (/^>$tag\s+/) ? eval {$s=$c} : 1;
 $c++;
}

print "$tag is at $s, and the file is $c. We have done ", int(($s/$c)*100000)/1000, " percent\n";
