#!/usr/bin/perl -w

# get DNA out of a genbank file

use Bio::SeqIO;
use strict;

my $usage=<<EOF;
$0 <list of genbankfiles>

EOF

die $usage unless ($ARGV[0]);

foreach my $file (@ARGV)
{
	my $of = $file;
	$of =~ s/\.gbk/.fasta/;
	$of =~ s/\.genbank/.fasta/;
	if ($of eq $file) {$of .= ".fasta"}
        my $sin=Bio::SeqIO->new(-file=>$file, -format=>'genbank');
	my $sout = Bio::SeqIO->new(-file=>">$of", -format=>'fasta');
	while (my $seq = $sin->next_seq()) {
		$sout->write_seq($seq);
	}
}

