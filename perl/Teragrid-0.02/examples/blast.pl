#!/usr/bin/perl -w

# A whole new blast program

use strict;
use lib '/homes/redwards/perl/share/perl/5.8.4/';
use TeraGrid::LSGW;
use Rob;
my $rob=new Rob;

my ($f, $n, $nf)=(undef, undef, 1);
while (@ARGV)
{
	my $t=shift;
	if ($t eq "-f") {$f=shift}
	if ($t eq "-n") {$n=shift}
	if ($t eq "-nf") {$nf=shift}
}

my $usage=<<EOF;
$0 
-f file name of fasta file
-n number of sequences to submit
-nf number of files to use (default = 1)

EOF

die $usage unless ($f && $n);

my $fa=$rob->read_fasta($f);

# randomize the sequences
my $keys=$rob->rand([keys %$fa]);
@$keys=splice(@$keys, 0, $n);

my @sequences;
my $inc=int($n/$nf)+1; # the +1 just means that there won't be 1 sequence
for (my $c=0; $c<=$#$keys; $c+=$inc)
{
	my $e=$c+$inc-1;
	my $content;
	foreach my $seq (@$keys[$c..$e])
	{
		next unless ($seq && $fa->{$seq});
		$content .= ">$seq\n". $fa->{$seq}. "\n";
	}       
	push @sequences, $content;
}       
print STDERR "Generated ", scalar(@sequences), " sequences to submit\n";



my $tg=new TeraGrid::LSGW(-verbose=>2);
my $c=0;
foreach my $sequence (@sequences)
{
	$c++;
	my $result=$tg->submit_blast(-sequence=>$sequence, -program=>"blastx");
	print STDERR "$c: $result\n";
}



