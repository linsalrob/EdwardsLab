# Schedule::SGE::Run

# POD docs

=head1 Schedule::SGE::Run

Submit jobs to the Sun Grid Engine (tm, probably), and check on the status of the engine. You should not use this method directly, rather you should use the Schedule::SGE method that inherits from this, then all the methods herein are available to you.

=head1 AUTHOR

 Rob Edwards (rob@salmonella.org)
 3/24/05

=cut

package Schedule::SGE::Run;
use strict;
use Exporter;

use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Schedule::SGE Exporter);
@EXPORT_OK = qw(command execute environment name project output_file error_file use_cwd notify mailto job_id);
our $VERSION = '0.01';

=head2 command()

Get or set the command that will be queued to by run

=cut

sub command {
 my ($self, @jobs)=@_;
 if ($jobs[0]) {
  $self->{'command'} = join " ", @jobs;
  $self->{'ran'}=0;
 }
 return $self->{'command'};
}


=head2 execute()

Run the command and return the queue id number.

The queue number is the location of the job in the queue, and you can cehck on that with queue_status

=cut

sub execute {
 my ($self) = @_;
 #return $self::Run->_run();
 return &_run($self);
}


=head2 environment()

Get and set the environment variables for Schedule::SGE. The variables that we need to care about are, and the default variables for my system are:
 	SGE_CELL		orionmulti
	SGE_EXECD_PORT		537
	SGE_QMASTER_PORT	536
	SGE_ROOT		/opt/sge

my $hashref=$sge->environment(\%vars);

=cut

sub environment {
 my ($self, $hash)=@_;
 my $return;
 foreach my $var (qw[SGE_CELL SGE_EXECD_PORT SGE_QMASTER_PORT SGE_ROOT]) {
  if ($hash->{$var}) {$ENV{$var}=$hash->{$var}}
  $return->{$var}=$ENV{$var};
 }
 return $return;
}
 
=head2 name()

Get or set the name of the job used by Schedule::SGE. 

=cut

sub name {
 my ($self, $val)=@_;
 if ($val) {
  unless ($val =~ /^[a-zA-Z]/) {
   print STDERR "Name must start with a letter. Name is now a$val\n";
   $val="a".$val;
  }
  if ($val =~ / /) {
   $val =~ s/ /_/g;
   print STDERR "Name can not have spaces in it. Name is now $val\n";
  }
  if (length($val) > 10) { 
   $val=substr($val, 0, 10);
   print STDERR "Name is truncated to 10 letters. Name is now $val\n";
  }
  $self->{'name'}=$val;
 }
 return $self->{'name'};
}

=head2 project()

Get or set the project used by Schedule::SGE. 

=cut

sub project {
 my ($self, $val)=@_;
 if ($val) {
  $val =~ s/ /_/g;
  $self->{'project'}=$val;
 }
 return $self->{'project'};
}

=head2 output_file()

Get or set the filename that will be used for the STDOUT

=cut

sub output_file {
 my ($self, $val)=@_;
 if ($val) {
  $self->{'output_file'}=$val;
 }
 return $self->{'output_file'};
}


=head2 error_file()

Get or set the filename that will be used for the STDERR

=cut

sub error_file {
 my ($self, $val)=@_;
 if ($val) {
  $self->{'error_file'}=$val;
 }
 return $self->{'error_file'};
}


=head2 use_cwd()

Boolean whether to set the cwd directory. NOTE: By default this is set to true, and you have to turn it off if you don't want it.

=cut

sub use_cwd {
 my ($self, $val)=@_;
 if (!defined $self->{'cwd'}) {$self->{'cwd'}=1}
 if (defined $val) {
  $self->{'cwd'}=$val;
 }
 return $self->{'cwd'};
}

=head2 notify()

Get or set whether the notify flag is set

=cut

sub notify {
 my ($self, $val)=@_;
 if (defined $val) {
  $self->{'notify'}=$val;
 }
 return $self->{'notify'};
}

=head2 mailto()

Email address to send the notify mail to

=cut

sub mailto {
 my ($self, $val)=@_;
 if ($val) {
  $self->{'mailto'}=1;
 }
 return $self->{'cwd'};
}


=head2 job_id()

The ID of the job that is submitted. This is only available after the command has begun, and is the ID of your job in the queue. Returns false if the job has not been executed or there was an error with the execution.

=cut

sub job_id {
 my ($self)=@_;
 return $self->{'job_id'};
}

=head2 _run()

An internal method to execute the command

=cut

sub _run {
 my ($self)=@_;
 return if ($self->{'ran'}); # we already run
 $self->{'ran'}=1;
 
 my $pipe = $self->executable('qsub');
 &_dieout('qsub') unless ($pipe);
 &_dieout('command') unless ($self->{'command'});
 
 my %tags=(
  'name' 	=> ' -N ',
  'project'	=> ' -P ',
  'mailto'	=> ' -M ',
  'output_file'	=> ' -o ',
  'error_file'	=> ' -e ',
 );

 foreach my $tag (keys %tags) {if ($self->{$tag}) {$pipe .= $tags{$tag}.$self->{$tag}}}

 if ($self->use_cwd) 		{$pipe .= " -cwd"}
 if ($self->notify)		{$pipe .= " -notify"}


 my $command=$self->{'command'};
 if ($self->{'verbose'}) {print STDERR "OPENING PIPE: $pipe\nSENDING JOB THROUGH PIPE: $command\n"}
 open(QSUB, "|$pipe > /tmp/$$.out 2>&1") || die "Can't open the pipe to submit jobs to";
 print QSUB $command, "\n";
 close QSUB;

 return 0 unless (-e "/tmp/$$.out"); 

 open (IN, "/tmp/$$.out") || die "Can't open /tmp/$$.out even though everything appeared to work fine";
 my $line=<IN>;
 close IN;
 $line =~ /Your job (\d+)/i;
 my $jobnumber=$1;
 if ($jobnumber) {
  unlink("/tmp/$$.out");
  $self->{'job_id'}=$jobnumber;
  return $jobnumber;
 }
 else {
  print STDERR "WARNING: No job number. This message was received:\n$line";
  return 0;
 }
}


=head2 _dieout()

Die nicely, with some kind of warning

=cut

sub _dieout {
 my ($self, $val)=@_;
 if ($val eq "command") {
  print STDERR <<EOF;

  You did not specify a command to run and so we are cowardly quitting.

EOF
 }
 elsif ($val eq "qsub") {
  print STDERR <<EOF;

  $0 could not find a $val executable. Please check to make sure that it is in your path, and you are running SGE.

EOF
 }
 else {
  print STDERR <<EOF;

  $0 died out for an unexplained reason. Sorry.

EOF
 }

 exit(-1);
}



1;
