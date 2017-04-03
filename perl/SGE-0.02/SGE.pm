# Schedule::SGE

# POD docs

=head1 Schedule::SGE

Interact with the Sun Grid Engine. This module locates the executables (qstat, qsub, etc), sets the level of verbosity, and other general commands related to the using the SGE.

=head1 ABSTRACT

Schedule::SGE is a suite of modules for interacting with the Sun Grid Engine. The base module Schedule::SGE handles locating the executables and making sure everything works fine. The three modules Schedule::SGE::Run, Schedule::SGE::Control, and Schedule::SGE::Status are for different interactions with the queues

=head1 AUTHOR

 Rob Edwards (rob@salmonella.org)
 3/24/05

=cut

package Schedule::SGE;
use strict;

use Schedule::SGE::Run qw/command execute environment name project output_file error_file use_cwd notify mailto job_id/;
use Schedule::SGE::Status qw/user status brief_job_stats all_jobs/;
use Schedule::SGE::Control qw/qdel/;

BEGIN {
 # set some default environment variables. These are default for me, but we'll provide a method to handle them too
 !$ENV{'SGE_CELL'} 		&& ($ENV{'SGE_CELL'}="orionmulti");
 !$ENV{'SGE_EXECD_PORT'} 	&& ($ENV{'SGE_EXECD_PORT'}="537");
 !$ENV{'SGE_QMASTER_PORT'} 	&& ($ENV{'SGE_QMASTER_PORT'}="536");
 !$ENV{'SGE_ROOT'} 		&& ($ENV{'SGE_ROOT'}="/opt/sge");
 our $VERSION=0.01;
};

=head2 new()

Instantiate the object, and preload some data. For example

my $sge=Schedule::SGE->new(
 -project=>'redwards',
 -mailto=>'rob@salmonella.org',
 -executable=>{qsub=>'/usr/local/bin/qsub', qstat=>'/usr/local/bin/qstat'},
 -verbose=>1,
);

=cut

sub new {
 my ($class, %args)=@_;
 my $self = bless{},$class;
 
 # just start with verbosity = 0. This avoids an error in the >= checks
 $self->{'verbose'}=0;

 # load up if we know what we have
 foreach my $key (keys %args) {
  my $nodash=$key;
  $nodash =~ s/^\-//;
  if ($nodash eq "executable") {
   $self->executable($args{$key});
  }
  elsif ($nodash eq "use_cwd") {
   $self->use_cwd($args{$key});
  }
  elsif ($nodash eq "name") {$self->name($args{$key})}
  else {
   $self->{$nodash}=$args{$key};
  }
 }
 return $self;
}


=head2 verbose()

Increase the level of error messages that are printed out.

=cut

sub verbose {
 my ($self, $val)=@_;
 if (defined $val) {
  $self->{'verbose'}=$val;
 }
 return $self->{'verbose'};
}


=head2 executable()

Get or set the executables that we will use. This method takes upto two arguments. With no arguments we will try and guess the settings that we need, and if we fail we will die. With a single argument we will return that executable path/program, guess it if we don't know it, and then finally fail. With two arguments we will assume that the second is the location of the executable (incl. path) of the first.

We will also take a reference to a hash as the single argument. In this case, we will use the hash as locations of the executables.

e.g.s:

# using a hash to set all the executables at once (recommended as we don't have to guess anything)
my $exec={'qsub'=>'/usr/local/bin/qsub', 'qstat'=>'/usr/local/bin/qstat'}
$sge->exectuable($exec);
my $pid=$sge->job_id;


# guessing all the executables (not recommended)
$sge->exectuables();
my $pid=$sge->job_id;

# getting the value for qsub
my $qsubexec=$sge->executable('qsub');

# setting a single value for qsub only
my $qsubexec=$sge->executable('qsub', '/usr/local/bin/qsub');

At the moment we try and figure out locations for each of the following applications
qstat
qsub
qdel

=cut

sub executable {
 my ($self, $exec, $path)=@_;

 print STDERR "Trying to find $exec\n" if ($exec && $self->{'verbose'} >= 2); 
 # these are the executables that we care about:
 my @want=(qw[qsub qstat qdel]);

 # first, if it is a reference add it
 if (ref($exec) eq "HASH") {
  foreach my $e (@want) {
   if ($exec->{$e}) {$self->{'execute'}->{$e}=$exec->{$e}}
  }
  return $self->{'execute'};
 }

 # now if we have both arg/var
 elsif ($exec && $path) {
  $self->{'execute'}->{$exec}=$path;
  return $path;
 }

 # now if we only have arg
 elsif ($exec) {
  if ($self->{'execute'}->{$exec}) {return $self->{'execute'}->{$exec}}
 }

 # if we get here we need to guess the locations 
 foreach my $e (@want) {
  my $guess=`which $e`; chomp($guess);
  if ($self->{'verbose'} >=2 ) {
   print STDERR "Looking for $e and found $guess\n";
  }
  
  if ($guess) {
   $self->{'execute'}->{$e}=$guess;
  } else {
   &_dieout("trying to get $e");
  }
 }

 if ($exec) {return $self->{'execute'}->{$exec}}
 else {return $self->{'execute'}}
}
  


1;
