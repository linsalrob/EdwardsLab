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


# we merge if the number and title is the same.
my %same = (
	"CS600" =>   ["METHDS BIOINFORMATICS MED", "BIOMI600"],
	"BIOL568" => ["BIOINFORMATICS", "BIOMI568"],
	"CS609" =>   ["COMPUTATIONAL GENOMICS", "BIOMI609"],
	"BIOL596" => ["ECOLOGICAL METAGENOMICS", "BIOL562"],
	"CS596" =>   ["PARALLEL COMPUTING", "CS505"],
	"CS596" =>   ["COMPUTATIONAL GENOMICS", "BIOMI609"]
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
	if ($same{$a[3]} && $same{$a[3]}[0] eq $a[4]) {$a[3] = $same{$a[3]}[1]}

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
		my $av = sprintf("%.2f", $rob->mean($count->{$c}));
		print join("\t", $c, $n, $av, $gotcredit->{$c}), "\n";
		delete $count->{$c};
	} else {
		print join("\t", $c, $gotcredit->{$c}, 0, $gotcredit->{$c}), "\n";
	}
}


foreach my $c (keys %$count) {
	my $n = $#{$count->{$c}} + 1;
	my $av = sprintf("%.2f", $rob->mean($count->{$c}));
	print join("\t", $c, $n, $av, 0), "\n";
}
