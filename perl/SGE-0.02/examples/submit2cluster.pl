#!/usr/bin/perl -w

# submit a job using my new method

use strict;
use Schedule::SGE;

my $usage=<<EOF;

$0 <options> <command>

OPTIONS: 
 -N name    (default = first word of command)
 -P project (default = redwards)
 
 -n notify  		(notify is not set unless -n is included. The default email address is redwards\@salmonella.org, but is unaffected unless you set notify).
 -e email address
 
EOF

my ($name, $project, $command, $notify, $email)=('', 'redwards', '', 0, "redwards\@salmonella.org");
while (@ARGV) {
 my $test=shift @ARGV;
 if ($test eq "-N") {$name=shift @ARGV}
 elsif ($test eq "-P") {$project=shift @ARGV}
 elsif ($test eq "-n") {$notify=1}
 elsif ($test eq "-e") {$email=shift @ARGV}
 else {$command .= " ". $test}
}

unless ($name) { 
 $command =~ m/^(\S+)/; $name=$1;
}
 

my $sge=Schedule::SGE->new(
 -project   	=> $project,
 -executable 	=> {qsub=>'/opt/sge/bin/lx26-x86/qsub', qstat=>'/opt/sge/bin/lx26-x86/qstat'},
 -name		=> $name,
 -verbose   	=> 1,
 -notify 	=> $notify,
 -mailto	=> $email,
);

$sge->command($command);

my $pid=$sge->execute();
exit(0);
my $stats;
until ($stats->[0]) {
 sleep(5);
 $stats=$sge->brief_job_stats($pid);
}

print "Stats\n", join "\t", @$stats, "\n";
 

