#!/usr/bin/perl -w

# Join two lists from two files. Anything with the same keys will be merged. This assumes there are two columns separated by tabs, and the first column is the key and the second the value.

use strict;

my $header; my @files; my $filetitle; my $zero = ""; my $skip=0;
while (@ARGV) {
	my $t=shift @ARGV;
	if ($t eq "-h") {$header=1}
	elsif ($t eq "-t") {$filetitle=1}
	elsif ($t eq "-z") {$zero=0}
	elsif ($t eq "-s") {$skip=1}
	elsif (-e $t) {push @files, $t}
	else {
		print STDERR "Don't understand $t\n";
	}
}

unless (scalar(@files) > 1) {
 print STDERR <<EOF;

$0 <files>

Join two or more lists from two files. Anything with the same keys will be merged. 
This assumes there are two separated by tabs, and the first column is the key and 
the rest of the columns the values

-h files include header row (first column is used from file 1)
-t use the file names as titles in the output
-z use 0 instead of null for non-existent values
-s skip lines that start #

EOF

die;
}

my $data; my %allkeys; my $headers; 
my %datapoints;
my $firstcolheader;
foreach my $f (@files) {
	open(IN, $f) || die "Can't open $f";
	$datapoints{$f}=0;
	while (<IN>) {
		chomp;
		if ($skip && index($_, "#") == 0) {next}
		my @a=split /\t/;
		my $key = shift @a;
		if ($header && !(defined $firstcolheader)) {$firstcolheader = $key}
		if ($header && !(defined $headers->{$f})) {$headers->{$f}=\@a}
		else {$data->{$f}->{$key} = \@a; $allkeys{$key}=1}
		($#a > $datapoints{$f}) ? $datapoints{$f} = $#a : 1;
	}
	close IN;
}

if ($header) {
    print $firstcolheader;
    map {print join("\t", "", @{$headers->{$_}})} @files;
    print "\n";
}

my @keys = sort {$a cmp $b} keys %allkeys;
if ($filetitle) {print join("\t", "", @files), "\n"}
foreach my $k (@keys) {
	print $k;
	foreach my $f (@files) {
		#(defined $data->{$f}->{$k}) ? (print join("\t", "", @{$data->{$f}->{$k}})) : print "\t" x scalar(@files);
		#(defined $data->{$f}->{$k}) ? (print join("\t", "", @{$data->{$f}->{$k}})) : print "\t" x $datapoints{$f};
		if (!defined $data->{$f}->{$k}) {
			$data->{$f}->{$k} = [];
			$#{$data->{$f}->{$k}}=$datapoints{$f};
		}
		map {(!defined $data->{$f}->{$k}->[$_]) ? $data->{$f}->{$k}->[$_]=$zero :1} (0 .. $#{$data->{$f}->{$k}});
		print join("\t", "", @{$data->{$f}->{$k}});
	}
	
	print "\n";
}


