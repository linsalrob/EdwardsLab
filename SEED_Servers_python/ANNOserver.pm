package ANNOserver;

#
# This is a SAS Component
#

=head1 Annotation Support Server Object


This file contains the functions and utilities used by the Annotation Support Server
(B<anno_server.cgi>). The various methods listed in the sections below represent
function calls direct to the server. These all have a signature similar to the
following.

    my $results = $annoObject->function_name($args);

where C<$annoObject> is an object created by this module, 
C<$args> is a parameter structure, and C<function_name> is the Annotation Support
Server function name. The output $results is a scalar, generally a hash
reference, but sometimes a string or a list reference.

=head2 Constructor

Use

    my $annoObject = ANNOserver->new();

to create a new annotation support server function object. The function object
is used to invoke the L</Primary Methods> listed below. See L<SAPserver> for
more information on how to create this object and the options available.

=cut

use LWP::UserAgent;
use Data::Dumper;
use YAML;

use base qw(ClientThing);

use strict;

sub new
{
    my($class, @options) = @_;
    my %options = ClientThing::FixOptions(@options);
    $options{url} = ClientThing::ComputeURL($options{url}, 'anno_server.cgi', 'anno');

#    my $self = {
#	server_url => $server_url,
#	ua => LWP::UserAgent->new(),
#    };
#    $self->{ua}->timeout(20 * 60);

#    return bless $self, $class;

    my $self = $class->SUPER::new(ANNO => %options);
    return $self;
}

#
# Doc for stuff in ANNO.pm is actually here.
#

=head1 Primary Methods

=head2 Functions

=head3 metabolic_reconstruction

    my $results = $annoObject->metabolic_reconstruction({
                                -roles => { [$role1, $id1],
                                            [$role2, $id2],
                                            ... });

This method will find for each subsystem, the subsystem variant that contains a
maximal subset of the roles in an incoming list, and output the ID of the
variant and a list of the roles in it.

=over 4

=item parameters

The single parameter is a reference to a hash containing the following
key fields.

=over 8

=item -roles

Reference to a list of 2-tuples, each consisting of the functional role followed
by an arbitrary ID of the caller's choosing (e.g., a gene name, a
sequence-project gene ID, a protein ID or whatever).

=back

For backward compatibility, instead of a hash reference you may specify a
simple reference to a list of 2-tuples.

=item RETURN

Returns a list of tuples, each containing a variant name (subsystem name,
colon, variant code), a role ID, and optionally a caller-provided ID associated
with the role. The ability to specify arbitrary IDs to be associated with the
roles is normally used to associate arbitrary gene IDs with the roles they are
believed to implement.

=back

=head3 find_special_proteins

    my $proteinList =   $annoObject->find_special_proteins({
                                -contigs => [[$contigID1, $contigNote1, $contigDNA1],
                                             [$contigID2, $contigNote2, $contigDNA2],
                                             ...],
                                -is_init => [$codon1a, $codon1b, ...],
                                -is_alt => [$codon2a, $codon2b, ...],
                                -is_term => [$codon3a, $codon3b, ...],
                                -comment => $commentString,
                                -templates => [[$protID1, $protNote1, $protSeq1],
                                              [$protID2, $protNote2, $protSeq1],
                                              ...]
                        });

This method searches for special proteins in a list of contigs. The method is
specifically designed to find selenoproteins and pyrrolysoproteins, but custom
protein templates can be specified to allow searching for any type of protein
family.

=over 4

=item parameter

The parameter is a reference to a hash with the following permissible keys.

=over 8

=item -contigs

Reference to a list of contigs. Each contig is represented by a 3-tuple
consisting of a contig ID, a comment, and a DNA string.

=item -is_init (optional)

Reference to a list of DNA codons to be used as start codons. The default is
C<ATG> and C<GTG>.

=item -is_alt (optional)

Reference to a list of DNA codons to be used as alternative start codons. These are
employed if there are no results from the main start codons. The default is
C<TTG>.

=item -is_term (optional)

Reference to a list of DNA codons to be used as stop codons. The default is
C<TAA>, C<TAG>, and C<TGA>.

=item -templates (optional)

Description of the type of special protein being sought. If C<pyrrolysoprotein>, then
the method will search for pyrrolysines. If C<selenoprotein>, then the method will
search for selenoproteins. Otherwise, should be a reference to a list of 3-tuples
containing templates for the proteins in the desired family. Each 3-tuple must
consist of an ID, a functional role description, and a protein sequence. The default
is C<selenoprotein>.

=item -comment (optional)

A string that will be inserted as a comment in each element of the output list. The
default is either C<pyrrolysoprotein> or C<selenoprotein>, depending on the template
specification.

=back

=item RETURN

Returns a reference to a list of hashes. Each hash contains the following keys.

=over 8

=item location

A location string describing the contig, start, and end location of the protein
found.

=item sequence

The protein sequence found.

=item reference_id

ID of the relevant template protein sequence.

=item reference_def

Functional role of the relevant template protein sequence.

=item comment

Comment from the input parameters.

=back

=back


=head3 assign_function_to_prot

    my $resultHandle = $annoObject->assign_function_to_prot($args)

For each incoming protein sequence, attempt to assign a function.
There are two ways functions can get assigned.  The first is based on
kmers, and these are normally viewed as the most reliable (at least
they give a consistent vocabulary!).  If no kmer match is made,
you can optionally try to make an assignment based on similarity
computations.

The attempt is made using kmer-technology.  A pass through the sequence
will locate "signature kmers", and scores will be computed.
The scores are based on the number of nonoverlapping hits, the number of
overlapping hits, and the difference in counts between hits against the
most probable function's kmer-set and the next most probable function's 
kmer set.  Basically, we compute all matching kmers.  Then, we split them
into sets based on the predictions each would make (each kmer, in effect,
predicts a single function).  One threshhold (the B<scoreThreshold>) is
the difference between total number of overlapping hits for the "best function"
versus the total number for the "next best".  B<hitTheshold> is the number
of overlapping hits required for the "best function".  Similarly,
B<seqHitThreshold> is the minimum number of non-overlapping hits.

Now, to add complexity, these thresholds are based on counting "1" for
each matched Kmer.  That is, the B<scoreThreshold> is normally thought of
as a difference in the number of occurrences.  However, you may wish to
"normalize" this threshold by dividing the count by the length of the sequence.
This then gives scores between 0 and 1, rather than between 0 and the length
of the sequence (-K if you wish to be pedantic).

=over 4

=item args

Reference to a hash containing the parameters. The allowable parameter fields
are as follows.

=over 8

=item -input

Either (1) an open input handle to a file containing the proteins in FASTA format,
or (2) A reference to a list of sequence data entries. Each entry is a triple of strings
[sequence-id, comment, protein-sequence-data].

=item -kmer

Specify the kmer size to use for analysis (valid sizes are 7 - 12).

=item -assignToAll

If TRUE, then if the standard matching algorithm fails to assign a protein,
a similarity-based assignment algorithm will be used instead.    

=item -scoreThreshold N

Require a Kmer score of at least N for a Kmer match to succeed.

=item -hitThreshold N

Require at least N (possibly overlapping) Kmer hits for a Kmer match to succeed.

=item -seqHitThreshold N

Require at least N (non-overlapping) Kmer hits for a Kmer match to succeed.

=item -normalizeScores 0|1

Normalize the scores to the size of the protein.

=item -detailed 0|1

If true, return a detailed accounting of the kmers hit for each protein.

=back

=item RETURN

Returns a Result Handle. Call C<get_next> on the result handle to get
back a data item. Each item sent back by C<get_next> is a 7-tuple
containing the results. Each tuple is of the form

    [ sequence-id, assigned-function, genome-set-name, score, non-overlapping hit-count, overlapping hit-count, detailed-hits]

where detailed-hits is undef unless the -detailed option was used.

If details were requested, the detailed-hit list is a list of tuples,
one for each kmer hit. These tuples have the form

    [ offset, oligo, functional-role, genome-set-name]

=back

=cut

sub assign_function_to_prot
{
    my($self, $args) = _handle_args(@_);

    my $wq;

    # my $params = [blast => $blast, min_hits => $min_hits, assign_to_all => ($assignToAll ? 1 : 0)];

    my $input = delete $args->{-input};
    
    my $params = [ map { $_ => $args->{$_} } keys %$args ];

    if (ref($input) eq 'ARRAY')
    {
	$wq = SequenceListWorkQueue->new($input);
    }
    else
    {
	$wq = FastaWorkQueue->new($input);
    }

    my $req_bytes = 1_000_000;
    #my $req_bytes = $blast ? 1000 : 1_000_000;

    return ResultHandler->new($self, $wq, $self->{server_url}, 'assign_function_to_prot', \&id_seq_pair_bundler,
			      #\&tab_delimited_output_parser,
			      \&YAML::Load,
			      $params, $req_bytes);
}

=head3 call_genes

    my $result = $annoObject->call_genes($args);

Call the protein-encoding genes for the specified DNA sequences. 

=over 4

=item args

Reference to a hash containing the parameters. The allowable parameter fields
are as follows.

=over 8

=item -input

The DNA sequences to be analyzed. This may take one of two forms: (1) a file handle
that is open for reading from a file of DNA sequences in FASTA format, or
(2) a reference to a list of DNA data entries. Each entry is a triple of strings
[sequence-id, comment, dna-sequence-data].

=item -trainingLocations (optional)

Reference to a hash mapping gene IDs to location strings. The location strings in
this case are of the form I<contig>C<_>I<start>C<_>I<end>. The locations indicated
should be coding regions in the incoming sequences to be analyzed (or in the
training contigs if a I<-trainingContigs> parameter is specified). If this
parameter is omitted, then the default GLIMMER training algorithm will be used.

=item -trainingContigs (optional)

The contigs in which the I<-trainingLocations> can be found. This may take one of
two forms: (1) a file handle that is open for reading from a file of DNA
sequences in FASTA format, or (2) a reference to a list of contig data entries.
Each entry is a triple of strings [contig-id, comment, contig-sequence-data].

=item -minContigLen (optional)

Shortest-length contig considered to be valid. This is used to prevent attempting
to call genes in contigs too short to have any complete ones. The default is C<2000>.

=item -geneticCode (optional)

The numeric code for the mapping from DNA to amino acids. The default is C<11>,
which is the standard mapping and should be used in almost all cases. A complete
list of mapping codes can be found at
L<http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>.

=back

=item RETURN

Returns a 2-tuple consisting of 1) a string containing what would normally
be the contents of an entire FASTA file for all the proteins found
followed by 2) a reference to a list of genes found. Each gene found will be
represented by a 4-tuple containing an ID for the gene, the ID of the contig
containing it, the starting offset, and the ending offset.

=back

=cut

use constant CALL_GENES_FILE_PARMS => { -input => 'id_seq', -trainingContigs => 'train_seq' };

sub call_genes
{
    my($self, $args) = _handle_args(@_);

    my $params = {
                  -geneticCode => ($args->{-geneticCode} || 11),
                  -minContigLen => ($args->{-minContigLen} || 2000),
                  -verbose => ($args->{-verbose} || 0),
                  };
    
    for my $fileKey (keys %{CALL_GENES_FILE_PARMS()}) {
        my $fileData = $args->{$fileKey};
        if (defined $fileData) {
            my $input = [];
            if (ref($fileData) eq 'ARRAY') {
                $input = $fileData;
            } else {
                my $fh;
                if (ref($fileData)) {
                    $fh = $fileData;
                } else {
                    my $fasta_file = $fileData;
                    open($fh, "<", $fasta_file);
                }
                while (my($id, $seqp, $com) = FastaWorkQueue::read_fasta_record($fh)) {
                    push(@$input, "$id,$$seqp");
                }
                close($fh);
            }
            $params->{CALL_GENES_FILE_PARMS->{$fileKey}} = $input;
        }
    }
    my $trainingData = $args->{-trainingLocations};
    if (ref $trainingData eq 'HASH') {
        $params->{training_set} = [ map { $_, $trainingData->{$_} } keys %$trainingData ];
    }
    my $parmList = [function => "call_genes", %$params];
    return $self->run_query_form($parmList);
}

=head3 find_rnas

    my $document = $annoObject->find_rnas($args)

Call the RNAs for the specified DNA sequences.

=over 4

=item args

Reference to a hash containing the parameters. The allowable parameter fields
are as follows.

=over 8

=item -input

Open input handle to a file containing the DNA sequences in FASTA format.

=item -genus

Common name of the genus for this DNA.

=item -species

Common name of the species for this DNA.

=item -domain

Domain of this DNA. The default is C<Bacteria>.

=back

=item RETURN

Returns a 2-tuple consisting of 1) a string containing what would normally
be the contents of an entire FASTA file for all the RNA genes found
followed by 2) a reference to a list of RNA genes found. Each gene found will be
represented by a 5-tuple containing an ID for the gene, the ID of the contig
containing it, the starting offset, the ending offset, and the
name of the RNA found.

=back

=cut

sub find_rnas
{
    my($self, $args) = _handle_args(@_);

    my $input = delete $args->{-input};
    
    my $params = [ map { $_ => $args->{$_} } keys %$args ];

    if (ref($input) ne 'ARRAY')
    {
	my $fh;
	if (ref($input))
	{
	    $fh = $input;
	}
	else
	{
	    my $fasta_file = $input;
	    open($fh, "<", $fasta_file);
	}
	$input = [];
	while (my($id, $seqp, $com) = FastaWorkQueue::read_fasta_record($fh))
	{
	    push(@$input, "$id,$$seqp");
	}
	close($fh);
    }

    return $self->run_query_form([function => "find_rnas",
				  @$params,
				  id_seq => $input]);
}

=head3 assign_functions_to_dna

    my $result = $annoObject->assign_functions_to_dna($args)

Analyze DNA sequences and output regions that probably belong to FIGfams.
The selected regions will be high-probability candidates for protein
encoding sequences.

=over 4

=item args

Reference to a hash containing the parameters. The allowable parameter fields
are as follows.

=over 8

=item -input

The sequences to be analyzed. This may take one of two forms:

1. An file handle that is open for reading from a file of DNA sequences in FASTA format, or
    
2. A reference to a list of sequence data entries. Each entry is a triple of strings
[sequence-id, comment, dna-sequence-data].

=item -kmer

Specify the kmer size to use for analysis  (valid sizes are 7 - 12).

=item -minHits

A number from 1 to 10, indicating the minimum number of matches required to
consider a protein as a candidate for assignment to a FIGfam. A higher value
indicates a more reliable matching algorithm; the default is C<3>.

=item -maxGap

When looking for a match, if two sequence elements match and are closer than
this distance, then they will be considered part of a single match. Otherwise,
the match will be split. The default is C<600>.

=back

=item RETURN

Returns a Result Handle. Call C<get_next> on the result handle to get back a data
item. Each item sent back by the result handle is a 2-tuple containing the
incoming protein sequence and a reference to a list of hit regions. Each hit
region is a 5-tuple consisting of the number of matches to the function, the start
location, the stop location, the proposed function, and the name of the
Genome Set from which the gene is likely to have originated.

=back

=cut

sub assign_functions_to_dna
{
    my($self, $args) = _handle_args(@_);
    my $wq;
    
    my $input = delete $args->{-input};
    
    my $params = [ map { $_ => $args->{$_} } keys %$args ];

    if (ref($input) eq 'ARRAY')
    {
	$wq = SequenceListWorkQueue->new($input);
    }
    else
    {
	$wq = FastaWorkQueue->new($input);
    }

    my $req_bytes = 500_000;

    return ResultHandler->new($self, $wq, $self->{server_url}, 'assign_functions_to_DNA',
			      \&id_seq_pair_bundler,
			      \&tab_delimited_output_parser, $params, $req_bytes);
}

###### Utility Methods ######

sub run_query
{
    my($self, $function, @args ) = @_;
    my $form = [function  => $function,
		args => YAML::Dump(\@args),
		];
    return $self->run_query_form($form);
}

sub run_query_form
{
    my($self, $form, $raw) = @_;

    my $res = $self->_send_request($form);
    #my $res = $self->{ua}->post($self->{server_url}, $form);
    
    if ($res)
    {
	my $content = $res;
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
	die "run_query_form: error encountered";
    }
}

sub id_seq_pair_bundler
{
    my($item) = @_;
    my($id, $seq) = @$item[0,2];
    return "id_seq", join(",", $id, (ref($seq) eq 'SCALAR' ? $$seq : $seq));
}

sub tab_delimited_output_parser
{
    my($line) = @_;
    chomp $line;
    my @cols = split(/\t/, $line);
    return \@cols;
}


sub tab_delimited_dna_data_output_parser
{
    my($line) = @_;
    chomp $line;
    my ($id, $idbe, $fam) = split(/\t/, $line);
    my ($beg, $end) = $idbe =~ /_(\d+)_(\d+)$/;
    return [$id, $beg, $end, $fam];
}


#
# Turn an argument list into a $self ref and an argument hash.
# Code lifted from ClientThing.
#
sub _handle_args
{
    my $self = shift;
    my $args = $_[0];
    if (defined $args)
    {
        if (scalar @_ gt 1)
	{
            # Here we have multiple arguments. We check the first one for a
            # leading hyphen.
            if ($args =~ /^-/) {
                # This means we have hash-form parameters.
                my %args = @_;
                $args = \%args;
            } else {
                # This means we have list-form parameters.
                my @args = @_;
                $args = \@args;
            }
        } else {
            # Here we have a single argument. If it's a scalar, we convert it
            # to a singleton list.
            if (! ref $args) {
                $args = [$args];
            }
        }
    }
    return($self, $args);
}




package ResultHandler;
use strict;
use Data::Dumper;

sub new
{
    my($class, $server_obj, $work_queue, $server_url, $function, $input_bundler, $output_parser, $form_vars, $req_bytes) = @_;

    my $self = {
	server_obj => $server_obj,
	work_queue => $work_queue,
	server_url => $server_url,
	function => $function,
	input_bundler => $input_bundler,
	output_parser => $output_parser,
	ua => LWP::UserAgent->new(),
	cur_result => undef,
	form_vars => $form_vars ? $form_vars : [],
	req_bytes => ($req_bytes ? $req_bytes : 16000),
    };
    $self->{ua}->timeout(20 * 60);
    return bless $self, $class;
}

sub get_next
{
    my($self) = @_;

    my $res =  $self->get_next_from_result();
    # print "gnfr returns: " , Dumper($res);

    if ($res)
    {
	return $res;
    }
    else
    {
	
	while (my @inp = $self->{work_queue}->get_next_n_bytes($self->{req_bytes}))
	{
	    my $form = [@{$self->{form_vars}}];
	    push(@$form, function => $self->{function},
			 map { &{$self->{input_bundler}}($_) } @inp);
	    #print "Invoke " .Dumper($form);

#	    my $res = $self->{ua}->post($self->{server_url}, $form);
	    my $res = $self->{server_obj}->_send_request($form);
	    if (defined($res))
	    {
		eval { 
		    $self->{cur_result} = [YAML::Load($res)];
		};
		if ($@)
		{
		    die "Query returned unparsable content ($@): " . $res->content;
		}
		#print "res: " . Dumper($self->{cur_result});
		my $oneres =  $self->get_next_from_result();
		if ($oneres)
		{
		    return $oneres;
		}
	    }
	    else
	    {
		die "error on post";
	    }
	}
	return;
    }
}

sub get_next_from_result
{
    my($self) = @_;
    my $l = $self->{cur_result};
    if ($l and @$l)
    {
	return shift(@$l);
    }
    else
    {
	delete $self->{cur_result};
	return undef;
    }
}

package SequenceWorkQueue;
use strict;

sub new
{
    my($class) = @_;

    my $self = {};
    
    return bless $self, $class;
}

sub get_next_n
{
    my($self, $n) = @_;
    my @out;
    
    for (my $i = 0;$i < $n; $i++)
    {
	my($id, $com, $seqp) = $self->get_next();
	if (defined($id))
	{
	    push(@out, [$id, $com, $seqp]);
	}
	else
	{
	    last;
	}
    }
    return @out;
}

sub get_next_n_bytes
{
    my($self, $n) = @_;
    my @out;

    my $size = 0;
    while ($size < $n)
    {
	my($id, $com, $seqp) = $self->get_next();
	if (defined($id))
	{
	    push(@out, [$id, $com, $seqp]);
	    $size += (ref($seqp) eq 'SCALAR') ? length($$seqp) : length($seqp);
	}
	else
	{
	    last;
	}
    }
    return @out;
}

package FastaWorkQueue;
use strict;
use base 'SequenceWorkQueue';
use FileHandle;

sub new
{
    my($class, $input) = @_;

    my $fh;
    if (ref($input))
    {
	$fh = $input;
    }
    else
    {
	$fh = new FileHandle("<$input");
    }

    my $self = $class->SUPER::new();

    $self->{fh} = $fh;

    return bless $self, $class;
}

sub get_next
{
    my($self) = @_;

    my($id, $seqp, $com) = read_fasta_record($self->{fh});
    return defined($id) ? ($id, $com, $seqp) : ();
}

sub read_fasta_record {
    my ($file_handle) = @_;
    my ($old_end_of_record, $fasta_record, @lines, $head, $sequence, $seq_id, $comment, @parsed_fasta_record);

    if (not defined($file_handle))  { $file_handle = \*STDIN; }

    $old_end_of_record = $/;
    $/ = "\n>";

    if (defined($fasta_record = <$file_handle>)) {
        chomp $fasta_record;
        @lines  =  split( /\n/, $fasta_record );
        $head   =  shift @lines;
        $head   =~ s/^>?//;
        $head   =~ m/^(\S+)/;
        $seq_id = $1;
        if ($head  =~ m/^\S+\s+(.*)$/)  { $comment = $1; } else { $comment = ""; }
        $sequence  =  join( "", @lines );
        @parsed_fasta_record = ( $seq_id, \$sequence, $comment );
    } else {
        @parsed_fasta_record = ();
    }

    $/ = $old_end_of_record;

    return @parsed_fasta_record;
}

package SequenceListWorkQueue;
use strict;
use base 'SequenceWorkQueue';

sub new
{
    my($class, $input) = @_;

    my $fh;
    if (ref($input) ne 'ARRAY')
    {
	die "SequenceWorkQueue requires a list as input";
    }

    my $self = $class->SUPER::new();

    $self->{list} = $input;

    return bless $self, $class;
}

sub get_next
{
    my($self) = @_;

    my $top = shift @{$self->{list}};

    return defined($top) ? @$top : ();
}


1;

