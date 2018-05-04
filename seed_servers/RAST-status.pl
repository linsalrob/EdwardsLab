#!/usr/bin/env perl

use strict;
use RASTserver;
$ENV{SAS_SERVER}="SEED";
use Term::ReadKey;

print "Please enter your RAST username:  ";
my $user = ReadLine(0);
chomp $user;

print "Please enter your RAST password:  ";
ReadMode 2;
my $password = ReadLine(0);
chomp $password;
ReadMode 1;
print "\n";


my $rast=new RASTserver($user, $password);
unless (defined $rast) {die "Can't connect ot the rast server"}





die "$0 <list of jobs>" unless (defined $ARGV[0]);
my $stat = $rast->status_of_RAST_job({-job => \@ARGV});

foreach my $job (sort {$a <=> $b} keys %$stat) {
	print join("\t", $job, $stat->{$job}->{'status'}), "\n";
}
