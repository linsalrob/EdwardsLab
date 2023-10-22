use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;

my $file = shift || die "file to process";
open(IN, $file) || die "$! : $file";
my @a = split /\t/, <IN>;
chomp(@a);
my $idx = shift @a;

my @first;
my @second;
map {if (m/chip/) {push @first, $_} elsif ($_) {push @second, $_}} @a;
my $c=0;
foreach my $l (sort {$a cmp $b} @first) {
	$c++;
	print "$l\t$c\n";
}
foreach my $l (sort {$a cmp $b} @second) {
	$c++;
	print "$l\t$c\n";
}


