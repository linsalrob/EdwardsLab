# TeraGrid::LSGW::Control

# POD docs

=head1 TeraGrid::LSW::Control

Access control and confirmation for the TeraGrid Life Sciences Gateway. This module handles instantiating the connection and ensuring that the username and password is correct. This should be called before submitting jobs to the queues, checking their status, or retrieving results.

=head1 AUTHOR

 Rob Edwards (RobE@theFIG.info)
 06/08/06

=head1 LICENSE

This software is released under the same terms as the SEED Toolkit. You can redistribute it and/or modify it under the terms of the SEED Toolkit Public License. 

You should have received a copy of the SEED Toolkit Public License along with this program; if not write to the University of Chicago at info@ci.uchicago.edu or the Fellowship for Interpretation of Genomes at veronika@thefig.info or download a copy from http://www.theseed.org/LICENSE.TXT.

=cut

package TeraGrid::LSGW::Control;
use strict;
use Exporter;
use Term::ReadKey; # for the password stuff
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
use HTTP::Cookies;

use vars qw(@ISA @EXPORT_OK);
@ISA = qw(TeraGrid::LSGW Exporter);
@EXPORT_OK = qw(connection login user password cookie_file);
our $VERSION = '0.02';

=head1 connection()

Return the connection to the TeraGrid gateway, and if one is not present, instantiate one

=cut

sub connection {
	my ($self)=@_;
	unless ($self->{'ua'}) {$self->{'ua'}=$self->login()}
	return $self->{'ua'};
}

=head1 login()

Login to the TeraGrid site using the username and password

Returns the useragent with the active connection if the login is successful, else returns undef

use: my $ua=$lsgw->login(username, password);

Note that if username and password are not provided they will be requested.

=cut

sub login {
	my ($self, $username, $password)=@_;
	$username=$self->user($username);
	$password=$self->password($password);

	$self->{'ua'} = LWP::UserAgent->new;
	$self->{'ua'}->agent("TG_LSGW_Rob/0.1 ");
	my $cf=$self->cookie_file();
	$self->{'ua'}->cookie_jar(HTTP::Cookies->new(file => $cf, autosave => 1));

	my $req = POST 'http://lsgw.mcs.anl.gov/applications_and_tools',
	   [
	   	'edit[name]'=> $username,
	   	'edit[pass]' => $password,
	   	'edit[form_id]' => 'user_login_block',
	   	'op' => 'Log in',
	   ];

	my $res = $self->{'ua'}->request($req);
	if ($res->status_line =~ /^500/) {die "Sorry, the webserver returned ", $res->status_line, "\n"}
	if ($self->{verbose}) {print STDERR "Logging in returned ", $res->status_line, "\n"}
	if ($res->status_line eq "302 Found") 
	{
		return $self->{'ua'}
	}
	else 
	{
		print STDERR "Sorry could not login with the information you gave\n";
		exit(-1);
	}
}


=head1 user()

Set the username. If a username is not provided then we'll request one

=cut

sub user {
	my ($self, $user)=@_;
	if ($user) {$self->{'user'}=$user}
	unless (defined $self->{'user'})
	{
		print "\nWe need a user authorized to submit jobs\nPlease enter the username:  ";
		$self->{'user'} = ReadLine 0;
		chomp $self->{'user'};
	}
	return $self->{'user'};
}

=head1 password()

Set the password. If a password is not provided then we'll request one

=cut

sub password {
	my ($self, $pass)=@_;
	if ($pass) {$self->{'password'}=$pass}
	unless (defined $self->{'password'})
	{
		print "Please enter the password :  ";
		ReadMode 2;
		$self->{'password'} = ReadLine 0;
		chomp $self->{'password'};
		ReadMode 1;
		print "\n";
	}
	return $self->{'password'};
}

=head1 cookie_file

Get and set the file to store the cookies in. At the moment the default is /tmp/$$.cookies

=cut

sub cookie_file {
	my ($self, $cf)=@_;
	if ($cf) {$self->{'cookie_file'}=$cf}
	unless (defined $self->{'cookie_file'}) {$self->{'cookie_file'}="/tmp/$$.cookies"}
	return $self->{'cookie_file'};
}


1;
