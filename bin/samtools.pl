#!/bin/env perl

# This is completely copied from Dave Tang's blog: https://davetang.org/muse/2014/03/06/understanding-bam-flags/
# please use the version there (I probably edited and broke this)
# and please give him credit and citations


use strict;
use warnings;

while (<>) {
	chomp;
	my $flag = $_;
	die "Please enter a numerical value\n" if $flag =~ /\D+/;

	my @output;

	if ($flag & 0x1){
		push @output, "multiple segments";
	}
	if ($flag & 0x2){
		push @output, "properly aligned";
	}
	if ($flag & 0x4){
		push @output, "segment unmapped";
	}
	if ($flag & 0x8){
		push @output, "next segment unmapped";
	}
	if ($flag & 0x10){
		push @output, "reverse complement";
	}
	if ($flag & 0x20){
		push @output, "next segment reverse complemented";
	}
	if ($flag & 0x40){
		push @output, "first segment";
	}
	if ($flag & 0x80){
		push @output,  "last segment";
	}
	if ($flag & 0x100){
		push @output,  "secondary alignment";
	}
	if ($flag & 0x200){
		push @output,  "not passing QC";
	}
	if ($flag & 0x400){
		push @output,  "PCR or optical duplicate";
	}
	if ($flag & 0x800){
		push @output,  "supplementary alignment";
	}
	print "$flag: ", join(" | ", @output),  "\n";
}
exit(0);
