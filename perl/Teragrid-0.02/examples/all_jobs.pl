#!/usr/bin/perl -w

# just print out a list of all jobs

use strict;
use lib '/homes/redwards/perl/share/perl/5.8.4/';
use TeraGrid::LSGW;


my $tg=new TeraGrid::LSGW(-verbose=>2);
my $aj=$tg->jobs();
print STDERR "There are ", scalar(keys %$aj), " jobs\n";
my $jl=$tg->job_list;
print join("\n", @$jl), "\n";
