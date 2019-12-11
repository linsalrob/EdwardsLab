use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;
my $rob = new Rob;

my $counter = 0;

my $usage = <<EOF;
$0 <id.map file> <output file> <list of fasta files>
id.map and output file must not exist (will not be automatically overwritten)
list of fasta files can be long
EOF

my $idf = shift || die $usage;
my $out = shift || die $usage;

if (-e $idf) {die "$idf exists.\n$usage"}
if (-e $out) {die "$out exists.\n$usage"}

open(IDF, ">$idf") || die "Can't write $idf";
open(OUT, ">$out") || die "Can't write $out";

foreach my $f (@ARGV) {
	my $fa = $rob->read_fasta($f);
	foreach my $seqid (keys %$fa) {
		$counter++;
		print OUT ">$counter\n", $fa->{$seqid}, "\n";
		print IDF "$f\t$seqid\t$counter\n";
	}
}

close(IDF);
close(OUT);
