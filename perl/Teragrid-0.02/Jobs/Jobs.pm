# TeraGrid::LSGW::Jobs

# POD docs

=head1 TeraGrid::LSGW::Jobs

Retrieve status, input, output, and results for jobs.

=head1 AUTHOR

 Rob Edwards (RobE@theFIG.info)
 06/08/06

=head1 LICENSE

This software is released under the same terms as the SEED Toolkit. You can redistribute it and/or modify it under the terms of the SEED Toolkit Public License. 

You should have received a copy of the SEED Toolkit Public License along with this program; if not write to the University of Chicago at info@ci.uchicago.edu or the Fellowship for Interpretation of Genomes at veronika@thefig.info or download a copy from http://www.theseed.org/LICENSE.TXT.

=cut

package TeraGrid::LSGW::Jobs;
use strict;
use Exporter;
use HTTP::Request::Common qw(GET POST);

use vars qw(@ISA @EXPORT_OK);
@ISA = qw(TeraGrid::LSGW Exporter);
@EXPORT_OK = qw(jobs input output results job_list);
our $VERSION = '0.02';


=head1 job_list()

Retrieve a list of jobs from the grid. This returns an array of jobs in the order in which they were submitted, with the most recent submissions first.

my $list=$lsgw->job_list();

=cut

sub job_list {
	my ($self)=@_;
	unless ($self->{'all_jobs'}) {$self->jobs()}
	return $self->{'all_jobs'};
}


=head1 jobs()

Retrieve a list of jobs, and their status, from the grid.

my $jobs=$lsgw->jobs();

$jobs is a reference to a hash. The keys are the jobs numbers and the values are a reference to an array of Name, Options, Input, Output, Results, Status

=cut

sub jobs {
	my ($self)=@_;
	my $cx=$self->connection();
	my $req = GET "http://lsgw.mcs.anl.gov/jobs";
	my $res=$cx->request($req);
	if ($res->status_line =~ /^500/) {die "Sorry, we got ", $res->status_line, "\n"}

	print STDERR $res->as_string, "\n";
	$res->as_string =~ /Job Summary Table.*?\<table(.*?)<\/table/is;
	my $results=$1;
	$results =~ s/\n//g;
	my @rows;
	while ($results =~ s/<tr.*?>(.*?)<\/tr.*?>//i) {push @rows, $1}
	foreach my $row (@rows)
	{
		my @cols;
		while ($row =~ s/<td.*?>(.*?)<\/td.*?>//i) {push @cols, $1}
		my $jid=shift @cols;
		push @{$self->{'all_jobs'}}, $jid;
		$self->{'job'}->{$jid}=\@cols;
	}
	return $self->{'job'};
}
	
=head1 input()

Retrieve the input to a job. This will return the data that was sent to the job (specifically the sequences) from the website

my $input = $lsgw->input($job);

=cut

sub input {
	my ($self, $job)=@_;
	unless ($self->{'job'}) {$self->jobs()}
	unless ($self->{'job'}->{$job})
	{
		print STDERR "No input is available for $job\n";
		return;
	}
	my $link=$self->{'job'}->{$job}->[2];
	$link =~ /(http:.*?)["'>]/;
	my $url=$1;
	my $cx=$self->connection();
	my $req = GET $url;
	return $cx->request($req)->as_string;
}

	
=head1 output()

Retrieve the output for a job. This will return the output from the teragrid PBS server

my $output = $lsgw->output($job);

=cut

sub output {
	my ($self, $job)=@_;
	unless ($self->{'job'}) {$self->jobs()}
	unless ($self->{'job'}->{$job})
	{
		print STDERR "No output is available for $job\n";
		return;
	}
	my $link=$self->{'job'}->{$job}->[3];
	$link =~ /(http:.*?)["'>]/;
	my $url=$1;
	unless ($url)
	{
		print STDERR "The output for $job does not have a url. Has it finished yet?\n";
		return;
	}
	my $cx=$self->connection();
	my $req = GET $url;
	return $cx->request($req)->as_string;
}

	
=head1 results()

Retrieve the results to a job. This will return the actual results of the job (i.e. the blast output)

my $results = $lsgw->results($job);

=cut

sub results {
	my ($self, $job)=@_;
	unless ($self->{'job'}) {$self->jobs()}
	unless ($self->{'job'}->{$job})
	{
		print STDERR "No results are available for $job\n";
		return;
	}
	my $link=$self->{'job'}->{$job}->[4];
	$link =~ /(http:.*?)["'>]/;
	my $url=$1;
	unless ($url)
	{
		print STDERR "The results for $job do not have a url. Has it finished yet?\n";
		return;
	}

	my $cx=$self->connection();
	my $req = GET $url;
	return $cx->request($req)->as_string;
}

1;
