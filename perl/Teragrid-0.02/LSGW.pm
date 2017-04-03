# TeraGrid::LSGW

# POD docs

=head1 TeraGrid::LSGW

Interact with the TeraGrid Life Sciences Gateway Interface using webservices

=head1 ABSTRACT

This is a set of perl modules for interacting with the Life Sciences Gateway interface to the TeraGrid. Initially written as an interface to the BLAST application it will be extended as needed to cover interpro too.

=head1 AUTHOR

Rob Edwards (RobE@theFIG.info)
06/08/06

=head1 LICENSE

This software is released under the same terms as the SEED Toolkit. You can redistribute it and/or modify it under the terms of the SEED Toolkit Public License. 

You should have received a copy of the SEED Toolkit Public License along with this program; if not write to the University of Chicago at info@ci.uchicago.edu or the Fellowship for Interpretation of Genomes at veronika@thefig.info or download a copy from http://www.theseed.org/LICENSE.TXT.

=cut


package TeraGrid::LSGW;
use strict;
use HTTP::Request::Common qw(POST);

use TeraGrid::LSGW::Control qw/connection login user password cookie_file/;
use TeraGrid::LSGW::Jobs qw/jobs input output results job_list/;

BEGIN {
 # set some default environment variables. These are default for me, but we'll provide a method to handle them too
 our $VERSION=0.02;
};

=head2 new()

Instantiate the object, and preload some data. For example

my $lsgw=TeraGrid::LSGW->new(
 -verbose=>1,
 -username=>$user,
 -password=>$password,
);

=cut

sub new {
 my ($class, %args)=@_;
 my $self = bless{},$class;
 
 # just start with verbosity = 0. This avoids an error in the >= checks
 $self->{'verbose'}=0;

 # load up if we know what we have
 foreach my $key (keys %args) {
  if ($key eq "-password") {$self->password($args{$key})}
  if ($key eq "-username") {$self->user($args{$key})}
  if ($key eq "-verbose") {$self->verbose($args{$key})}
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

  
=head2 submit_blast()

Submit a sequence to be compared to a library using BLAST.

Arguments: 
	
	-sequence => sequence string in fasta format joined by newline (reqd)
	-program  => BLAST program (reqd)
	-database => database to comapre to (default = Seed_NR_2006-04-27)
	-options  => BLAST options (default= -m 8 -e 0.1)

my $success=$lsgw->submit_blast(-sequence=>">fasta\nGATCAGCATGCCATG\n", -database=>"Seed_NR_2006-04-27");

=cut

sub submit_blast {
	my ($self, %args)=@_;
	my ($seq, $prog, $db, $opt)=(undef, undef, "Seed_NR_2006-04-27", "-m 8 -e 0.1");
	foreach my $k (keys %args)
	{
		if ($k eq "-sequence") {$seq=$args{$k}}
		elsif ($k eq "-program") {$prog=$args{$k}}
		elsif ($k eq "-database") {$db=$args{$k}}
		elsif ($k eq "-options") {$opt=$args{$k}}
	}

	unless ($seq) {print STDERR "Can't run blast, no sequence defined\n"; return}
	unless ($prog) {print STDERR "Can't run blast, no blast program defined\n"; return}

	my $req = POST 'http://lsgw.mcs.anl.gov/cgi-bin/run.app',
		Content_Type => 'form-data',
		Content =>
		[
			'submit' => 'BLAST',
			'user' => 7,
			'database' => $db,
			'blast-options' => $opt,
			'blastprogram' => 'blastx',
			'tool' => 'blast',
			'sequence' => $seq,
		];

	my $cx=$self->connection();
	my $res=$cx->request($req);
	if ($res->status_line =~ /^500/) {die "Sorry, we got ", $res->status_line, "\nHTML is\n", $res->as_string, "\n"}
	print STDERR "Submitting BLAST returned ", $res->status_line, "\n" if ($self->verbose());
	#print $res->as_string;
	return 1 if ($res->status_line eq "200 OK");
}
	



1;
