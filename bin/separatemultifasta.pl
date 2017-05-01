#!/usr/bin/perl -w

# rewritten with BioPerl

use strict;
use Bio::SeqIO;

my $file=shift || die "$0 <fasta file>";
my $dir=$file.".files";
if (-e $dir) {die "$dir already exists. Not overwriting"}
else {mkdir $dir, 0755}

my %seen;
my $sio=Bio::SeqIO->new(-file=>$file, -format=>"fasta");
while (my $seq=$sio->next_seq)
{
 my $id=$seq->id;
 $id =~ s/\s+/_/g;
 while ($seen{$id}) {print "$id already written, "; $id.="1"; print " now trying $id\n"}
 my $fout=Bio::SeqIO->new(-file=>">$dir/$id.fasta", -format=>"fasta");
 $seen{$id}=1;
 $fout->write_seq($seq);
}

