# Schedule::SGE::Status

# POD docs

=head1 Schedule::SGE::Status

Check on the status of the Schedule::SGE queues. You should not use this method directly, rather you should use the Schedule::SGE method that inherits from this, then all the methods herein are available to you.

=head1 AUTHOR

 Rob Edwards (rob@salmonella.org)
 3/24/05

=cut

package Schedule::SGE::Status;
use strict;
use Exporter;

use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Schedule::SGE Exporter);
@EXPORT_OK = qw(user status brief_job_stats all_jobs);
our $VERSION = '0.01';

=head2 user()

Set the user of the processes. If not defined will be guess by whoami

=cut

sub user {
 my ($self, $user)=@_;
 $self->{'user'} = $user if ($user);
 unless ($self->{'user'}) {
  $self->{'user'}=`whoami`;
  chomp($self->{'user'});
 }
 return $self->{'user'};
}
 
=head2 status()

Get the queue status. This will return a hash where each key is the name of a node, and each value is a reference to an array. The array has the following components:

0. Queue type (one of B(atch), I(nteractive), C(heckpointing), P(arallel), T(ransfer) or combinations thereof or N(one))
1. Processors used
2. Load average
3. State

=cut

sub status {
 my ($self)=@_;
 my $qstat=$self->executable('qstat');
 my @status = `$qstat -f`;
 my $node;
 my $allstats;
 while (@status) {
  my $line=shift @status;
  chomp($line);
  next if ($line =~ /^\-/);
  next if ($line =~ /queuename/);
  next if (!$line || $line =~ /^\s*$/);
  if ($line =~ /^\S+/ && $line !~ /^\#/) {
   # it is a node status
   my @pieces = split /\s+/, $line;
   $node=shift @pieces;
   my $type=shift @pieces;
   if ($type =~ /[^BICPTN]/) {
    die "Received type of $type from \n$line\n\tbut it should not contain anything other than B,I,C,P,T or N";
   }
   my $procs=shift @pieces;
   if ($procs !~ m#\d+/\d+#) {
    die "Received procs of $procs from \n$line\n\tbut it should not contain anything other than \\d+/\\d+";
   }
   my ($load_avg, $arch, $states)= @pieces;
   if ($self->verbose && $states) {print STDERR "Node $node has state $states\n"}
   unless (defined $states) {$states=''}
   $allstats->{$node}=[$type, $procs, $load_avg, $states];
  }
  elsif ($line =~ m#^\s+(\d+)\s+(\d+\.\d+)\s+(\S+.*?)\s+\S+\s+(\S+)\s+(\d+/\d+/\d+\s+\d+\:\d+\:\d+)\s+\S+#) {
   # it is a job
   # something like 
   # 1441 0.56000 testing123 rob          r     03/25/2005 11:59:12     1
   my ($pid, $load, $name, $user, $date)=($1, $2, $3, $4, $5);
   $self->{'job'}->{$pid}=[$node, $pid, $load, $name, $user, $date];
  }
  elsif ($line =~ /^\#/ || $line =~ /PENDING/) {
   # at the end of the list there are some pending jobs
   while (@status) {
    my $pend=shift(@status);
    next if ($pend =~ /PENDING/ || $pend =~ /^\#/);
    $pend =~ s/^\s+//; $pend =~ s/\s+$//;
    my @pieces=split /\s+/, $pend;
    next unless (scalar @pieces > 5);
    my ($user, $status, $date, $time, $processes)=splice(@pieces, -5, 5);
    my ($pid, $load)=splice(@pieces, 0, 2);
    $date .= " ".$time;
    my $name = join " ", @pieces;
    $self->{'job'}->{$pid}=['pending', $pid, $load, $name, $user, $date];
   }
  }
  else {
   print STDERR "We don't know how to parse |$line|\n";
  }
 }


 return $allstats;
}


=head2 brief_job_stats()

Get some brief statistics about a job. This method will return a reference to an array with the following statistics:
Node the job is/was running on
Process ID of the job
Load
Name of process
Username
Date and time of submission

for example, my $stats=$sge->brief_job_stats($job);

=cut

sub brief_job_stats {
 my ($self, $job)=@_;
 return [] if (!$job);
 return $self->{'job'}->{$job} if ($self->{'job'}->{$job});
 
 # if we get this far, we should run status and quickly get the status
 $self->status();
 if ($self->{'job'}->{$job}) {
  return $self->{'job'}->{$job};
 }
 else {
  return [];
 }
}


=head2 all_jobs()

Returns an array of all jobs that were found in the queues.

=cut

sub all_jobs {
 my ($self)=@_;
 unless ($self->{'job'}) {$self->status()}
 return keys %{$self->{'job'}};
}


1;
