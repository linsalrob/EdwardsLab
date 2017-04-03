#!/usr/bin/perl -w

# Test the jobs interface

use strict;
use lib '/homes/redwards/perl/share/perl/5.8.4/';
use TeraGrid::LSGW;


my $job=shift || die "$0 <job number>\n";

my $tg=new TeraGrid::LSGW(-verbose=>2);

open(OUT, ">$job.input") || die "Can't write to $job.input";
print OUT join("\n", $tg->input($job)), "\n";
close OUT;

open(OUT, ">$job.output") || die "Can't write to $job.output";
print OUT join("\n", $tg->output($job)), "\n";
close OUT;

open(OUT, ">$job.results") || die "Can't write to $job.results";
print OUT join("\n", $tg->results($job)), "\n";
close OUT;

