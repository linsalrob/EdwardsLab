#!/usr/bin/perl
use strict;
use LWP::Simple;
use Data::Dumper;

my $range=shift || die "$0 <range> <number to get (default=10)> <output file>\nFor the range use something like this example AAVX02000001-AAVX02067420\nOutput file is optional";
my $n = shift;
my $outf = shift;
chomp $n;
chomp $outf;
unless ($n) {$n=10}

if ($n > 100) {die "Sorry. That number ($n) is too many, and you won't get any results. Please try a smaller number. It should probably be about 100 or less\n"}


my ($from, $to)=split /\-/, $range;
$from =~ /^(\D*0*)(\d+)$/; my ($fromid, $beg)=($1, $2);
$to =~ /^(\D*0*)(\d+)$/; my ($toid, $end)=($1, $2);

unless ($fromid eq $toid) {die "Couldn't get the base id's right. Had |$fromid| and |$toid| but these are not the same\n"}
my $id=$fromid;

if ($end < $beg) {($beg, $end)=($end, $beg); ($from, $to)=($to, $from)}

print "About to retrieve ", $end-$beg, " sequences.\n";
#print "Continue [Y/n]\n";
#my $in=<STDIN>;
#chomp($in);
#exit if (lc($in) =~ /^n/);

unless (defined $outf) {$outf = "${id}_sequences.gbk"}

open(OUT, ">$outf") || die "Can't write to $outf";
my $time=time;
my $url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&retmode=text&rettype=gb&id=';
my $urlid=$url;
my $added=0;
for (my $i=$beg; $i<=$end; $i++)
{
        my $val=$id.$i;
        $urlid .= "$val,";
        $added++;
        unless ($added % $n)
        {
                print STDERR "Getting sequences upto $added\n";
                $urlid =~ s/,$//;
                while (time-$time < 3) {sleep 1}
                $time=time;
		# print STDERR $urlid;
                print OUT get($urlid);
                $urlid=$url;
        }
}

print STDERR "Getting sequences upto $added\n";
$urlid =~ s/,$//;
while (time-$time < 3) {sleep 1}
# print STDERR $urlid;
print OUT get($urlid);

exit 0;

