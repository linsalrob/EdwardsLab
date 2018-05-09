#!/usr/bin/env perl

use strict;
use RASTserver;
$ENV{SAS_SERVER}="SEED";
use Term::ReadKey;

die "$0 [-u username] [-p password] <list of jobs>" unless (defined $ARGV[0]);

my @ids;
my $user; my $password;
while (@ARGV) {
	my $t=shift @ARGV;
	if ($t eq "-u") {$user = shift @ARGV}
	elsif ($t eq "-p") {$password = shift @ARGV}
	else {push @ids, $t}
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

my $stat = $rast->status_of_RAST_job({-job => \@ids});

foreach my $job (sort {$a <=> $b} keys %$stat) {
	print join("\t", $job, $stat->{$job}->{'status'}), "\n";
}
