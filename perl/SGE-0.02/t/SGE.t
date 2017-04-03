# -*-Perl-*-
# Test script for Schedule::SGE modules 
#
# Written by Rob Edwards, rob@salmonella.org
# 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
my $error;
BEGIN {
    $NUMTESTS = 22;
    use Test;
    plan tests => $NUMTESTS;
    $error = 0;
}

if( $error ==  1 ) {
    exit(0);
}

require Schedule::SGE;

ok(1);

# test Schedule::SGE
my $sge=new Schedule::SGE;
ok $sge;
my $verbose=$sge->verbose(1);
ok $verbose;
my $exec=$sge->executable;
ok $exec;
ok $exec->{'qdel'};
ok $exec->{'qsub'};
ok $exec->{'qstat'};

# test Schedule::SGE::Control
# ok $sge->qdel('noone'); # this is tricky, because I donot want to delete jobs that people have running.


# test Schedule::SGE::Run
my $run=$sge->command();
ok !($run); # run should not actually return anything!
my $environment=$sge->environment();
ok $environment;
ok $environment->{SGE_CELL};
ok $environment->{SGE_ROOT};
ok $environment->{SGE_QMASTER_PORT};
ok $environment->{SGE_EXECD_PORT};
ok !($sge->name());
ok !($sge->project());
ok !($sge->output_file());
ok !($sge->error_file());

ok $sge->use_cwd;
ok !($sge->notify());
ok $sge->mailto;


# test Schedule::SGE::Status
my  $user=$sge->user();
#print STDERR "User is $user\n";
ok $user;

#my $stat=$sge->status();
#ok $stat;
ok $sge->brief_job_stats();
#ok $sge->all_jobs();



