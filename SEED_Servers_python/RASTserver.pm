package RASTserver;

#
# This is a SAS Component
#

use LWP::UserAgent;
#use LWP::Debug qw(+trace +debug +conns);
use Data::Dumper;
use YAML;
use File::Basename;
use ClientThing;
use Cwd 'abs_path';

$YAML::CompressSeries = 0;

use strict;


=head1 RAST Submission Server Function Object

This file contains the functions and utilities used by the RAST Submission Server
. The various methods listed in the sections below represent
function calls direct to the server. These all have a signature similar to the
following.

    my $document = $rastObject->function_name($args);

where C<$rastObject> is an object created by this module,
C<$args> is a parameter structure, and C<function_name> is the RAST
Server function name. The output document is a structure, generally a hash
reference, but sometimes a string or a list reference.


=head2 Constructor

Use

    my $rastObject = RASTserver->new($username, $password);

    Where $username and $password are your RAST username and password. 


=cut




sub new
{
    my($class, $username, $password, $opts) = @_;

    my $server_url;
    if ($opts->{-test})
    {
	$server_url = "http://servers.nmpdr.org/rast-test/server.cgi";
    }
    elsif ($opts->{-server})
    {
	$server_url = $opts->{-server};
    }

    my $self = {
	server_url => ClientThing::ComputeURL($server_url, 'rast_submit_server.cgi', 'rast'),
	ua => LWP::UserAgent->new(),
	username => $username,
	password => $password,
    };
    $self->{ua}->timeout(20 * 60);

    return bless $self, $class;
}

=head1 Primary Methods

=head3 submit_RAST_job

	my $jobH = $rastObject->submit_RAST_job({

		-taxonomyID => The NCBI taxonomy id   | Domain => [Archaea | Bacteria | Virus]
            	-filetype => [Fasta|Genbank]
		if Fasta && ! TaxonomyID
			-geneticCode =>  [4 | 11]
	        -organismName  => Name of the organism
	        -keepGeneCalls =>  [0 | 1], 
		-geneCaller    => [RAST | Glimmer3], 
		-file => full path to the input file
}
);

This method submits a job to the RAST server. It returns a hash of:

	    {jobid} =id
	    {status} = [submitted|error]
	    {error_message} = message

=cut

		
sub submit_RAST_job
{
    my ($self, $params) = @_;
    my ($result, $error_string);
    
    my $user     = $self->{username};
    my $file     = $params->{-file};
    my $filetype = lc($params->{-filetype});
    my $tax_id   = $params->{-taxonomyID};
    my $domain   = lc($params->{-domain});
    my $organism = $params->{-organismName};
    my $keep     = $params->{-keepGeneCalls};
    my $genetic_code = $params->{-geneticCode};
    my $gene_caller  = lc($params->{-geneCaller});

    my $is_euk_user  = defined($FIG_Config::rast_euk_users)
	&& defined($FIG_Config::rast_euk_users->{$user})
	&& $FIG_Config::rast_euk_users->{$user};
    
    my %allowed_domains = map { $_ => 1 } (qw(Bacteria Archaea Viruses), ($is_euk_user ? q(Eukaryota) : ()));
    my %allowed_codes   = map { $_ => 1 } (4, 11, ($is_euk_user ?  (1) : ()));
    
    if ($filetype ne "fasta" && $filetype ne "genbank") {
	$error_string .= "-filetype not Fasta or Genbank\n";
    }
    
    if ($tax_id) {
	if ($tax_id =~ /\D/) {
		$error_string .= "-taxonomyID must be numeric\n";
        }
    } elsif ($filetype ne 'genbank' ) {
	if (0 == grep { lc(substr($_,0,1)) eq lc(substr($domain,0,1)) } (keys %allowed_domains)) {
	    my @tmp  = sort keys %allowed_domains;
	    $tmp[-1] = q(or ) . $tmp[-1];
	    my $tmp  = join(q(, ), @tmp);
	    $error_string .= "-domain must be $tmp\n";
    	}
	
	if ($organism !~ m/^[A-Z][a-z]*\s+\S+/o) {
	    $error_string .= "-organismName is mandatory for non-GenBank submissions\n";
	}
    }
    
    if (!$file) {
	$error_string .= "You must supply a file\n";
    } else {
	if (! open( F, "<$file" )) {
	    $error_string .= "Invalid file path\n";
	}
    }
    if (!defined($keep) || ($keep != 0 && $keep != 1)) {
	$error_string .= "-keepGeneCalls must be 0 or 1\n";
    }
    
    if ($filetype eq "fasta" && ! $tax_id) {
	if (!$genetic_code || (not $allowed_codes{$genetic_code})) {
	    my @tmp  = sort { $a <=> $b } (keys %allowed_codes);
	    $tmp[-1] = q( or ) . $tmp[-1];
	    my $tmp  = join(q(, ), @tmp);
	    $error_string .= "-geneticCode must be $tmp\n";
	}
    }
    if ($gene_caller ne "rast" && $gene_caller ne "glimmer3") {
	$error_string .= "-geneCaller must be RAST or Glimmer3\n";
    }
    
    
    if ($error_string) {
	$result->{status} = 'error';
	$result->{error_message} = $error_string;
	return ($result);
    }
    
    {
    	local $/;
	undef $/;
	$params->{-file} = <F>;
    }
    close(F);

    return $self->run_query('submit_RAST_job', $params);
}

sub get_contig_ids_in_project_from_entrez
{
    my($self, $params) = @_;

    return $self->run_query('get_contig_ids_in_project_from_entrez', $params);
}

sub get_contigs_from_entrez
{
    my($self, $params) = @_;

    return $self->run_query('get_contigs_from_entrez', $params);
}


=head3 status_of_RAST_job

	my $rastjobH = $rastObject->status_of_RAST_job({-job => \@jobs});

Where @jobs is a list of jobs

The return value is a hash keyed by Jobid of   

		    {status} = Job Stage
		    {error_msg} = message
		    {verbose-status} = RAST Metadata

=cut



sub status_of_RAST_job
{
    my($self, $params) = @_;

    return $self->run_query('status_of_RAST_job', $params);
}


=head3 retrieve_RAST_job

	my $result = $rastObject->retrieve_RAST_job(-job => $jobid, -format => $format)

where $jobid is the RAST id of the job and

$format is  one of:

		genbank 		(Genbank format)
		genbank_stripped 	(Genbank with EC numbers removed)
		embl 			(EMBL format)
		embl_stripped 		(EMBL with EC numbers stripped)
		gff3 			(GFF3 format)
		gff3_stripped 		(GFF3 with EC numbers stripped)
		rast_tarball 		(gzipped tar file of the entire job)

The return is a hash of 

		{status} = ok|error
		{file} = the downloaded file name
		{error_msg} = The error message
	
=cut
 

sub retrieve_RAST_job
{
    my($self, $params) = @_;

    my $filehandle = delete $params->{-filehandle};

    my $cb = sub {
	my($chunk) = @_;
	print $filehandle $chunk;
    };
    
    my $form = [function  => 'retrieve_RAST_job',
		args => YAML::Dump($params),
		username => $self->{username},
		password => $self->{password},
		];

    my $res = $self->{ua}->post($self->{server_url}, $form,
				':content_cb' => $cb,
			       );

    if ($res->is_success)
    {
	return { status => 'ok' };
    }
    else
    {
	return { status => 'error', error_msg => $res->status_line . $res->content };
    }
}	     




=head3 kill_RAST_job


	my $ret = rastObject->kill_RAST_job(-job => \@jobids);

where @jobids is an array of RAST job ids to kill.


Return is a hash keyed by Job ID of 

			{status} = ok|error
			{messages} = Messages

=cut



sub kill_RAST_job
{
    my($self, $params) = @_;

    return $self->run_query('kill_RAST_job', $params);
}



=head3 delete_RAST_job


	my $ret = rastObject->delete_RAST_job(-job => \@jobids);

where @jobids is an array of RAST job ids to kill.


Return is a hash keyed by Job ID of 

			{status} = ok|error
			{error_message} = Error Message

=cut

sub delete_RAST_job
{
    my($self, $params) = @_;

    return $self->run_query('delete_RAST_job', $params);
}


=head3 get_job_metadata


	my $ret = rastObject->get_job_metadata($jobid)

where $jobid is the RAST id  of a RAST job


Return is a hash of 

			{status} = ok|error
			{error_message} = Error Message
			{key} => {metdata}
=cut



sub get_job_metadata
{
    my($self, $params) = @_;
    return $self->run_query('get_job_metadata', $params);
}

sub run_query
{
    my($self, $function, @args ) = @_;
    my $form = [function  => $function,
		args => YAML::Dump(@args),
		username => $self->{username},
		password => $self->{password},
		];
    return $self->run_query_form($form);
}

sub run_query_form
{
    my($self, $form, $raw) = @_;

    my $res = $self->{ua}->post($self->{server_url}, $form);
    
    if ($res->is_success)
    {
	my $content = $res->content;
	if ($raw)
	{
	    return $content;
	}
	     
#	print "Got $content\n";
	my $ret;
	eval {
	    $ret = Load($content);
	};
	if ($@)
	{
	    die "Query returned unparsable content ($@): " . $content;
	}
	return $ret;
    }
    else
    {
	die "error on post " . $res->status_line . " " . $res->content;
    }
}

=head3 copy_to_RAST_dir

	my $result = $rastObject->copy_to_RAST_dir(-job  => $jobid, 
						   -from => $file_or_dir,
						   -to   => $to_dir);

where $jobid is the RAST id of the job,

      $file_or_dir is a file or directory to be copied into the RAST UserSpace, and

      $to_dir is optional.  If omitted the copy is to UserSpace.  If present, the
            copy is to UserSpace/$to_dir (intermediate directories get created, if 
            necessary)

The return is a hash of 

		{status} = ok|error
		{error_msg} = The error message
	
=cut

sub copy_to_RAST_dir
{
    my($self, $params) = @_;

    my $from = delete $params->{-from};
    my $to_name;
    my $src_fh;
    
    if (-f $from)
    {
	if (!open($src_fh, "<", $from))
	{
	    return { status => 'error',
			 error_msg => "Could not open $from: $!" };
	}
	$params->{-type} = 'file';
    }
    elsif (-d $from)
    {
	if (!open($src_fh, "-|", "tar", "-c", "-f", "-", $from))
	{
	    return { status => 'error',
			 error_msg => "Could not open tar pipe for $from: $!" };
	}
	$params->{-type} = 'tar';
    }
    else
    {
	return { status => 'error',
		 error_msg => "From not found: $from"
	       };
    }

    my $block_size = 50_000_000;
    
    $to_name = basename($from);
    $params->{-toName} = $to_name;

    my $buf;
    my $chunk_num = 0;
    my $size = 0;

#    my @stat = stat($src_fh);
#    my $file_size = $stat[7];
    
    while (my $n = read($src_fh, $buf, $block_size))
    {
	$params->{-chunkNum} = $chunk_num;
	
	my $form = [function  => 'copy_to_RAST_dir',
		    args => YAML::Dump($params),
		    username => $self->{username},
		    password => $self->{password},
		    file => $buf,
		    ];

	warn "Sending chunk " . $chunk_num + 1 . "\n";
	my $res = $self->{ua}->post($self->{server_url},
				    Content_Type => 'form-data',
				    Content => $form,
				   );

	if ($res->is_success)
	{
	    my $content = $res->content;
	    my $ret;
	    eval {
		$ret = Load($content);
	    };
	    if ($@)
	    {
		die "Query returned unparsable content ($@): " . $content;
	    }
	    warn "response " . Dumper($ret);
	    if ($ret->{status} ne 'ok')
	    {
		return $ret;
	    }
	}
	else
	{
	    return { status => 'error', error_msg => $res->status_line };
	}
	$size += $n;
	$chunk_num++;
    }
    $params->{-totalSize} = $size;
    
    my $form = [function  => 'copy_to_RAST_dir',
		args => YAML::Dump($params),
		username => $self->{username},
		password => $self->{password},
		];
    
    warn "send total size $size\n";
    my $res = $self->{ua}->post($self->{server_url},
				Content_Type => 'form-data',
				Content => $form,
			       );
    
    if ($res->is_success)
    {
	my $content = $res->content;
	my $ret;
	eval {
	    $ret = Load($content);
	};
	if ($@)
	{
	    die "Query returned unparsable content ($@): " . $content;
	}
	warn "response " . Dumper($ret);
	return $ret;
    }
    else
    {
	return { status => 'error', error_msg => $res->status_line };
    }
}	     

1;

