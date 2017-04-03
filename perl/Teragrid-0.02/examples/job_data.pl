#!/usr/bin/perl -w

# Test the jobs interface

use strict;
use lib '/homes/redwards/perl/share/perl/5.8.4/';
use TeraGrid::LSGW;


my $ajob="blastx-20060608-20153622";
#my $ajob="blastx-20060612-01534283";

my $tg=new TeraGrid::LSGW(-verbose=>2);

my $jobs=$tg->jobs();

foreach my $j (keys %$jobs)
{
	print join("\t", $j, @{$jobs->{$j}}), "\n";
}

print "For job $ajob\n";
print "INPUT\n======\n", $tg->input($ajob), "\nOUTPUT\n======\n", $tg->output($ajob), "\nRESULTS\n======\n", $tg->results($ajob), "\n";


