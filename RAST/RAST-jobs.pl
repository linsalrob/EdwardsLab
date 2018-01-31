#!/usr/bin/perl -w
#
#

use strict;
use RASTserver;
use Term::ReadKey;
use Data::Dumper;
$ENV{SAS_SERVER}="PUBSEED";

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

my $time = time; my $job = 0;
my @jobs = $rast->jobs();

foreach my $j (@jobs) {
	print Dumper($j);
	print STDERR $job++, " : ", ($time-time), " seconds\n";
}
