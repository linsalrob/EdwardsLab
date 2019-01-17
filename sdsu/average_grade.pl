#__perl__
#
#

use strict;
use Rob;
my $rob = new Rob;
my %grade = (
	"A" =>4.00,
	"A-" =>3.70,
	"B+" =>3.33,
	"B" =>3.00,
	"B-" =>2.70,
	"C+" =>2.30,
	"C" =>2.00,
	"C-" =>1.70,
	"D+" =>1.30,
	"D" =>1.00,
	"D-" =>0.70,
	"F" => 0,
);

my $f = shift || die "File to open";
open(IN, $f) || die "$f";
my $count;
my $credit;
my $gotcredit;

# ADDED RP to grade as 0
#ADDED W to grade as 0
#ADDED F to grade as 0
#ADDED IC to grade as 0
#ADDED WU to grade as 0
#ADDED AU to grade as 0
#
while (<IN>) {
	chomp;
	my @a=split /\t/;
	unless ($a[6]) {
		# print STDERR "No grade found in $_\n";
		# These are usually classes in the POS that people are currently taking
		next;
	}
	if ($a[6] eq "RP" || $a[6] eq "W" || $a[6] eq "IC" || $a[6] eq "WU" || $a[6] eq "AU") {
		# print STDERR "Skipped $_\n";
		next;
	}

	if ($a[6] eq "CR" || $a[6] eq "NC") {
		$credit->{"$a[3]\t$a[4]"}++;
		unless (defined $gotcredit->{"$a[3]\t$a[4]"}) {$gotcredit->{"$a[3]\t$a[4]"}=0}
		if ($a[6] eq "CR") {
			$gotcredit->{"$a[3]\t$a[4]"}++;
		}
		next;
	}

	unless (defined $grade{$a[6]}) {
		print STDERR "ADDED $a[6] to grade as 0\n";
		$grade{$a[6]} = 0;
	}
	push @{$count->{"$a[3]\t$a[4]"}}, $grade{$a[6]};
}

print("Class\tTitle\tNo.Students\tAv. grade\tCredit\n");

foreach my $c (keys %$credit) {
	if (defined $count->{$c}) {
		my $n = $credit->{$c} + length($count->{$c});
		print join("\t", $c, $n, $rob->mean($count->{$c}), $gotcredit->{$c}), "\n";
		delete $count->{$c};
	} else {
		print join("\t", $c, $gotcredit->{$c}, 0, $gotcredit->{$c}), "\n";
	}
}


foreach my $c (keys %$count) {
	my $n = $#{$count->{$c}} + 1;
	print join("\t", $c, $n, $rob->mean($count->{$c}), 0), "\n";
}
