use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;
my $rob = new Rob;

my $reads;
my %alltypes;
my %allsamples;
foreach my $f (@ARGV) {

	open(IN, $f) || die "$! : $f";
	while (<IN>) {
		chomp;
		s/.001.fast..gz//;# replace fastq or fasta
		s/.fast..gz//;
		s#//#/#g;
		my @a=split /\t/;
		my @b = split /\//, $a[0];
		$reads->{$b[1]}->{$b[0]}=$a[1];
		$allsamples{$b[1]}++;
		$alltypes{$b[0]}++;
	}
	close IN;
}

my @types = sort {$a cmp $b} keys %alltypes;
my @samples = sort {$a cmp $b} keys %allsamples;

print("Sample\t", join("\t", @types), "\n");
foreach my $s (@samples) {
	print $s;
	foreach my $t (@types) {
		print "\t";
		(defined $reads->{$s}->{$t}) ? print $reads->{$s}->{$t} : print "not found";
	}
	print "\n";
}


