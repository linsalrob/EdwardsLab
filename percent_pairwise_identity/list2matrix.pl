#!/usr/bin/env perl
# convert a file with [id1, id2, score] to an id x id matrix. 
# assumes the score is a % and sets id1==id2 to 100

use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;
my $rob = new Rob;

my %opts;
getopts('f:vo:', \%opts);
unless ($opts{f} && $opts{o}) {
	die <<EOF;
	$0
	-f tsv file to parse (required)
	-o matrix file to write (required)
	-v verbose output
EOF
}

my $score;
open(IN, $opts{f}) || die "$! : $opts{f}";
while (<IN>) {
	chomp;
	next if (/^id1/);
	s/>//g; # fix fasta headers
	my @a=split /\t/;
	$score->{$a[0]}->{$a[1]} = $a[2];
	$score->{$a[1]}->{$a[0]} = $a[2];
}

my @all = sort {$a cmp $b} keys %$score;

open(OUT, ">$opts{o}") || die "Can't write ot $opts{o}";
print OUT "". join("\t", @all), "\n";
foreach my $a (@all) {
	print OUT $a;
	foreach my $b (@all) {
		if ($a eq $b) {
			print OUT "\t100";
			next;
		}
		if (!defined $score->{$a}->{$b}) {
			print STDERR "No score for $a $b\n";
			print OUT "\t";
			next;
		}
		print OUT "\t", $score->{$a}->{$b};
	}
	print OUT "\n";
}
close OUT;





