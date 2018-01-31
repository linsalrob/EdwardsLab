#!/usr/bin/env perl 

use strict;
use Data::Dumper;
use RASTserver;
$ENV{SAS_SERVER}="SEED";
use Term::ReadKey;



die <<EOF  unless (defined $ARGV[0]); 
$0
[-d directory to output them into]
[-u user name]
[-p password]
[-f format]

<list of jobs>
EOF

my @jobs;
my $dir;
my $user; my $password;
my $format = 'rast_tarball';
while (@ARGV) {
	my $t=shift @ARGV;
	if ($t eq "-d") {$dir = shift}
	elsif ($t eq "-u") {$user = shift}
	elsif ($t eq "-p") {$password = shift}
	elsif ($t eq "-f") {$format = shift}
	else {push @jobs, $t}
}

if (!$user) {
	print "Please enter your RAST username:  ";
	$user = ReadLine(0);
	chomp $user;
}
if (!$password) {
	print "Please enter your RAST password:  ";
	ReadMode 2;
	$password = ReadLine(0);
	chomp $password;
	ReadMode 1;
	print "\n";
}

my $rast=new RASTserver($user, $password);
unless (defined $rast) {die "Can't connect ot the rast server"}


$dir && (! -e $dir) && (mkdir $dir, 0755);
(!$dir) && ($dir = ".");

my $stat = $rast->status_of_RAST_job({-job => \@jobs});
	#print Dumper($stat);

my @complete;
foreach my $job (sort {$a <=> $b} keys %$stat) {
	if ($stat->{$job}->{'status'} eq "complete") {push @complete, $job}
	else {print STDERR "$job is not complete. Ignored\n"}
}


my $ext = 'tar';
if ($format =~ /genbank/) {$ext = 'gbk'}
elsif ($format =~ /embl/) {$ext = 'embl'}
elsif ($format =~ /gff3/) {$ext = 'gff3'}

foreach my $job (@complete) {
	my $OUT;
	open($OUT,  ">$dir/$job.$ext") || die "Can't write to $dir/$job.$ext";
	my $f = $rast->retrieve_RAST_job({-job => $job, -format=>$format, -filehandle => $OUT});
	print "Retrieved $job with status ", $f->{status}, ". Writing to $dir/$job.$ext \n";
	close $OUT;
}	

