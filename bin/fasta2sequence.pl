=pod

Extract one or more sequences from a fasta file

=cut

use strict;
use Data::Dumper;
use Rob;
my $rob = new Rob;

unless ($#ARGV >= 1) {
	die "$0 <fasta fiile> [list of sequences to extract]\n";
}

my $faf = shift;
my $fa = $rob->read_fasta($faf);
foreach my $seq (@ARGV) {
	if (!defined $fa->{$seq}) {
		print STDERR "ERROR: $seq not found in $faf\n";
		next;
	}
	print ">$seq\n", $fa->{$seq}, "\n";
}




