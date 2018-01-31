#!/usr/bin/perl -w
#

use strict;
use Data::Dumper;
$ENV{SAS_SERVER}="PUBSEED";
print STDERR "SAS is $ENV{SAS_SERVER}\n";
use Term::ReadKey;
use RASTserver;

## Use RAST test, not regular RAST
#my $rast=new RASTserver("redwards", "couchman", {"-test"=>1});
# Now using regular RAST


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

print Dumper($rast->jobs());
