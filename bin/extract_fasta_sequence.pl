=pod

Extract one or more sequences from a fasta file, but not using much memory

=cut

use strict;

unless ($#ARGV >= 1) {
	die "$0 <fasta fiile> [list of sequences to extract]\n";
}

my $faf = shift;
my %want;
map {$want{$_}=1} @ARGV;

if ($faf =~ /.gz$/) {
	open(IN, "gunzip -c $faf|") || die "Can't open a pipe to $faf";
} else {
	open(IN, $faf) || die "$! $faf";
}
my $p=0;
while (<IN>) {
	if (index($_, ">") == 0) {
		$p = 0;
		if (/^>(\S+)/ && $want{$1}) {$p=1}
	}
	print if ($p);
}
close IN;



