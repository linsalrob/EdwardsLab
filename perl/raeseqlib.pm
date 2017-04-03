package raeseqlib;

######################################################################
#
# PLEASE READ THIS:
#
# This code came from Gary Olsen and is really called gjoseqlib. You should use that version
# This is a fork made by Rob so he could hack away at somethings and not offend Gary (or anyone else)
# all of this code is (c) Gary Olsen, and you should not use this module. If you don't know
# where to get an original version, ask Rob or try www.theseed.org/
#
#


#  A sequence entry is ( $id, $def, $seq )
#  A list of entries is a list of references
#
#  Efficient reading of an entire file of sequences:
#
#  @seq_entries = read_fasta( )                     # STDIN
#  @seq_entries = read_fasta( \*FILEHANDLE )
#  @seq_entries = read_fasta(  $filename )
#
#  Reading sequences one at a time to conserve memory.  Calls to different
#  files can be intermixed.
#
#  @entry = read_next_fasta_seq( \*FILEHANDLE )
# \@entry = read_next_fasta_seq( \*FILEHANDLE )
#  @entry = read_next_fasta_seq(  $filename )
# \@entry = read_next_fasta_seq(  $filename )
#  @entry = read_next_fasta_seq()                   # STDIN
# \@entry = read_next_fasta_seq()                   # STDIN
#
#  Legacy interface:
#  @seq_entries = read_fasta_seqs( \*FILEHANDLE )   # Original form
#
#  Reading clustal alignment.
#
#  @seq_entries = read_clustal( )                   # STDIN
#  @seq_entries = read_clustal( \*FILEHANDLE )
#  @seq_entries = read_clustal(  $filename )
#
#  Legacy interface:
#  @seq_entries = read_clustal_file(  $filename )
#
#  $seq_ind   = index_seq_list( @seq_entries );   # hash from ids to entries
#  @seq_entry = seq_entry_by_id( \%seq_index, $seq_id );
#  $seq_desc  = seq_desc_by_id(  \%seq_index, $seq_id );
#  $seq       = seq_data_by_id(  \%seq_index, $seq_id );
#
#  ( $id, $def ) = parse_fasta_title( $title )
#  ( $id, $def ) = split_fasta_title( $title )
#
#  Write a fasta format file from sequences.
#
#  print_alignment_as_fasta(                @seq_entry_list ); # STDOUT
#  print_alignment_as_fasta(               \@seq_entry_list ); # STDOUT
#  print_alignment_as_fasta( \*FILEHANDLE,  @seq_entry_list );
#  print_alignment_as_fasta( \*FILEHANDLE, \@seq_entry_list );
#  print_alignment_as_fasta(  $filename,    @seq_entry_list );
#  print_alignment_as_fasta(  $filename,   \@seq_entry_list );
#
#  Legacy interface:
#  print_seq_list_as_fasta( \*FILEHANDLE, @seq_entry_list );  # Original form
#
#  Interface that it really meant for internal use to write the next sequence
#  to an open file:
#
#  print_seq_as_fasta( \*FILEHANDLE, $id, $desc, $seq );
#  print_seq_as_fasta(               $id, $desc, $seq );
#  print_seq_as_fasta( \*FILEHANDLE, $id,        $seq );
#  print_seq_as_fasta(               $id,        $seq );
#
#  Write PHYLIP alignment.  Names might be altered to fit 10 character limit:
#
#  print_alignment_as_phylip(                @seq_entry_list ); # STDOUT
#  print_alignment_as_phylip(               \@seq_entry_list ); # STDOUT
#  print_alignment_as_phylip( \*FILEHANDLE,  @seq_entry_list );
#  print_alignment_as_phylip( \*FILEHANDLE, \@seq_entry_list );
#  print_alignment_as_phylip(  $filename,    @seq_entry_list );
#  print_alignment_as_phylip(  $filename,   \@seq_entry_list );
#
#  Write basic NEXUS alignment for PAUP.
#
#  print_alignment_as_nexus(               [ \%label_hash, ]  @seq_entry_list );
#  print_alignment_as_nexus(               [ \%label_hash, ] \@seq_entry_list );
#  print_alignment_as_nexus( \*FILEHANDLE, [ \%label_hash, ]  @seq_entry_list );
#  print_alignment_as_nexus( \*FILEHANDLE, [ \%label_hash, ] \@seq_entry_list );
#  print_alignment_as_nexus(  $filename,   [ \%label_hash, ]  @seq_entry_list );
#  print_alignment_as_nexus(  $filename,   [ \%label_hash, ] \@seq_entry_list );
#
#  print_gb_locus( \*FILEHANDLE, $locus, $def, $accession, $seq );
#
#  Remove extra columns of alignment gaps from an alignment:
#
#   @packed_seqs = pack_alignment(  @seqs )
#   @packed_seqs = pack_alignment( \@seqs )
#  \@packed_seqs = pack_alignment(  @seqs )
#  \@packed_seqs = pack_alignment( \@seqs )
#
#  Pack mask for an alignment (gap = 0x00, others are 0xFF)
#
#   $mask = alignment_gap_mask(  @seqs )
#   $mask = alignment_gap_mask( \@seqs )
#
#  Pack a sequence alignment according to a mask:
#
#   @packed = pack_alignment_by_mask( $mask,  @align )
#   @packed = pack_alignment_by_mask( $mask, \@align )
#  \@packed = pack_alignment_by_mask( $mask,  @align )
#  \@packed = pack_alignment_by_mask( $mask, \@align )
#
#  Expand sequence by a mask, adding indel at "\000" or '-' in mask:
#
#   $expanded = expand_sequence_by_mask( $seq, $mask )
#
#  Remove all alignment gaps from sequences (modify a copy):
#
#   @packed_seqs = pack_sequences(  @seqs )  #  Works for one sequence, too
#   @packed_seqs = pack_sequences( \@seqs )
#  \@packed_seqs = pack_sequences(  @seqs )
#  \@packed_seqs = pack_sequences( \@seqs )
#
# Basic sequence manipulation functions:
#
#  @entry  = subseq_DNA_entry( @seq_entry, $from, $to [, $fix_id] )
#  @entry  = subseq_RNA_entry( @seq_entry, $from, $to [, $fix_id] )
#  $DNAseq = DNA_subseq(  $seq, $from, $to )
#  $DNAseq = DNA_subseq( \$seq, $from, $to )
#  $RNAseq = RNA_subseq(  $seq, $from, $to )
#  $RNAseq = RNA_subseq( \$seq, $from, $to )
#  @entry  = complement_DNA_entry( @seq_entry [, $fix_id] )
#  @entry  = complement_RNA_entry( @seq_entry [, $fix_id] )
#  $DNAseq = complement_DNA_seq( $NA_seq )
#  $RNAseq = complement_RNA_seq( $NA_seq )
#  $DNAseq = to_DNA_seq( $NA_seq )
#  $RNAseq = to_RNA_seq( $NA_seq )
#  $seq    = pack_seq( $sequence )        # modifies a copy
#  $seq    = clean_ae_sequence( $seq )
#
#  $aa = translate_seq( $nt, $met_start )
#  $aa = translate_seq( $nt )
#  $aa = translate_codon( $triplet );
#
#  User-supplied genetic code.  The supplied code needs to be complete in
#  RNA and/or DNA, and upper and/or lower case.  The program guesses based
#  on lysine and phenylalanine codons.
#
#  $aa = translate_seq_with_user_code( $nt, $gen_code_hash, $met_start )
#  $aa = translate_seq_with_user_code( $nt, $gen_code_hash )
#
#  Locations (= oriented intervals) are ( id, start, end )
#  Intervals are ( id, left, right )
#
#  @intervals = read_intervals( \*FILEHANDLE )
#  @intervals = read_oriented_intervals( \*FILEHANDLE )
#  @intervals = standardize_intervals( @interval_refs ) # (id, left, right)
#  @joined    = join_intervals( @interval_refs )
#  @intervals = locations_2_intervals( @locations )
#  $interval  = locations_2_intervals( $location  )
#  @reversed  = reverse_intervals( @interval_refs )      # (id, end, start)
#
#  Convert GenBank locations to SEED locations
#
#  @seed_locs = gb_location_2_seed( $contig, @gb_locs )
#
#  Read quality scores from a fasta-like file:
#
#  @seq_entries = read_qual( )               #  STDIN
# \@seq_entries = read_qual( )               #  STDIN
#  @seq_entries = read_qual( \*FILEHANDLE )
# \@seq_entries = read_qual( \*FILEHANDLE )
#  @seq_entries = read_qual(  $filename )
# \@seq_entries = read_qual(  $filename )
#
#  Evaluate alignments:
#
#  $fraction_diff = fraction_nt_diff( $seq1, $seq2, \%options )
#  $fraction_diff = fraction_nt_diff( $seq1, $seq2 )
#  $fraction_diff = fraction_nt_diff( $seq1, $seq2, $gap_weight )
#  $fraction_diff = fraction_nt_diff( $seq1, $seq2, $gap_open, $gap_extend )
#
#  ( $npos, $nid, $ndif, $ngap, $nopen, $tgap, $topen ) = interpret_nt_align( $seq1, $seq2 )
#  ( $npos, $nid, $ndif, $ngap, $nopen, $tgap, $topen ) = interpret_aa_align( $seq1, $seq2 )
#
#  @sims = oligomer_similarity( $seq1, $seq2, \%opts )
#
#===============================================================================

use strict;
use Carp;
use Data::Dumper;

#  Exported global variables:

our @aa_1_letter_order;  # Alpha by 1 letter
our @aa_3_letter_order;  # PAM matrix order
our @aa_n_codon_order;  
our %genetic_code;
our %genetic_code_with_U;
our %amino_acid_codons_DNA;
our %amino_acid_codons_RNA;
our %n_codon_for_aa;
our %reverse_genetic_code_DNA;
our %reverse_genetic_code_RNA;
our %DNA_letter_can_be;
our %RNA_letter_can_be;
our %one_letter_to_three_letter_aa;
our %three_letter_to_one_letter_aa;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(
        read_fasta_seqs
        read_fasta
        read_next_fasta_seq
        read_clustal_file
        read_clustal
        parse_fasta_title
        split_fasta_title
        print_seq_list_as_fasta
        print_alignment_as_fasta
        print_alignment_as_phylip
        print_alignment_as_nexus
        print_seq_as_fasta
        print_gb_locus

        index_seq_list
        seq_entry_by_id
        seq_desc_by_id
        seq_data_by_id

        pack_alignment
        alignment_gap_mask
        pack_alignment_by_mask
        expand_sequence_by_mask
        pack_sequences

        subseq_DNA_entry
        subseq_RNA_entry
        DNA_subseq
        RNA_subseq
        complement_DNA_entry
        complement_RNA_entry
        complement_DNA_seq
        complement_RNA_seq
        to_DNA_seq
        to_RNA_seq
        pack_seq
        clean_ae_sequence

        translate_seq
        translate_codon
        translate_seq_with_user_code

        read_intervals
        standardize_intervals
        join_intervals
        locations_2_intervals
        read_oriented_intervals
        reverse_intervals

        gb_location_2_seed

        read_qual

        fraction_nt_diff
        interpret_nt_align
        interpret_aa_align
        oligomer_similarity
        );

our @EXPORT_OK = qw(
        @aa_1_letter_order
        @aa_3_letter_order
        @aa_n_codon_order
        %genetic_code
        %genetic_code_with_U
        %amino_acid_codons_DNA
        %amino_acid_codons_RNA
        %n_codon_for_aa
        %reverse_genetic_code_DNA
        %reverse_genetic_code_RNA
        %DNA_letter_can_be
        %RNA_letter_can_be
        %one_letter_to_three_letter_aa
        %three_letter_to_one_letter_aa
        );



=head1 new

 instantiate for OO perl

=cut

sub new {
	my ($class) = @_;
	my $self = {};
	bless $self, $class;
	return $self
}


#-----------------------------------------------------------------------------
#  Helper function for defining an input filehandle:
#     filehandle is passed through
#     string is taken as file name to be openend
#     undef or "" defaults to STDOUT
#
#    ( \*FH, $name, $close [, $file] ) = input_filehandle( $file );
#
#-----------------------------------------------------------------------------
sub input_filehandle
{
    my $file = shift;

    #  FILEHANDLE

    return ( $file, $file, 0 ) if ( ref( $file ) eq "GLOB" );

    #  Null string or undef

    return ( \*STDIN, "", 0 ) if ( ! defined( $file ) || ( $file eq "" ) );

    #  File name

    if ( ! ref( $file ) )
    {
        my $fh;
        if    ( -f $file                       ) { }
        elsif (    $file =~ /^>(.+)$/ && -f $1 ) { $file = $1 }
        else { die "Could not find input file '$file'\n" }
        open( $fh, "<$file" ) || die "Could not open '$file' for input\n";
        return ( $fh, $file, 1 );
    }

    #  Some other kind of reference; return the unused value

    return ( \*STDIN, undef, 0, $file );
}


#-----------------------------------------------------------------------------
#  Read fasta sequences from a filehandle (legacy interface; use read_fasta)
#  Save the contents in a list of refs to arrays:  (id, description, seq)
#
#     @seq_entries = read_fasta_seqs( \*FILEHANDLE )
#-----------------------------------------------------------------------------
sub read_fasta_seqs { read_fasta( @_ ) }


#-----------------------------------------------------------------------------
#  Read fasta sequences.  Save the contents in a list of refs to arrays:
#
#     $seq_entry = [ id, description, seq ]
#
#     @seq_entries = read_fasta( )               #  STDIN
#    \@seq_entries = read_fasta( )               #  STDIN
#     @seq_entries = read_fasta( \*FILEHANDLE )
#    \@seq_entries = read_fasta( \*FILEHANDLE )
#     @seq_entries = read_fasta(  $filename )
#    \@seq_entries = read_fasta(  $filename )
#  #  @seq_entries = read_fasta( "command |" )   #  open and read from pipe
#  # \@seq_entries = read_fasta( "command |" )   #  open and read from pipe
#     @seq_entries = read_fasta( \$string )      #  reference to file as string
#    \@seq_entries = read_fasta( \$string )      #  reference to file as string
#
#-----------------------------------------------------------------------------
sub read_fasta
{
    my @seqs;
    if ( $_[0] && ref $_[0] eq 'SCALAR' )
    {
        @seqs = map { $_->[2] =~ tr/ \n\r\t//d; $_ }
                map { /^(\S+)([ \t]+([^\n]*\S)?\s*)?\n(.+)$/s ? [ $1, $3 || '', $4 ] : () }
                split /^>\s*/m, ${$_[0]};
    }
    else
    {
        @seqs = map { $_->[2] =~ tr/ \n\r\t//d; $_ }
                map { /^(\S+)([ \t]+([^\n]*\S)?\s*)?\n(.+)$/s ? [ $1, $3 || '', $4 ] : () }
                split /^>\s*/m, slurp( @_ );
    }

    wantarray() ? @seqs : \@seqs;
}

#-----------------------------------------------------------------------------
#  A fast file reader:
#
#     $data = slurp( )               #  \*STDIN
#     $data = slurp( \*FILEHANDLE )  #  an open file handle
#     $data = slurp(  $filename )    #  a file name
#     $data = slurp( "<$filename" )  #  file with explicit direction
#   # $data = slurp( "$command |" )  #  open and read from pipe
#
#  Note:  It is faster to read lines by reading the file and splitting
#         than by reading the lines sequentially.  If space is not an
#         issue, this is the way to go.  If space is an issue, then lines
#         or records should be processed one-by-one (rather than loading
#         the whole input into a string or array).
#-----------------------------------------------------------------------------
sub slurp
{
    my ( $fh, $close );
    if ( ref $_[0] eq 'GLOB' )
    {
        $fh = shift;
    }
    elsif ( $_[0] && ! ref $_[0] )
    {
        my $file = shift;
        if    ( -f $file                       ) { $file = "<$file" }
        elsif (    $file =~ /^<(.*)$/ && -f $1 ) { }  # Explicit read
      # elsif (    $file =~ /\S\s*\|$/         ) { }  # Read from a pipe
        else                                     { return undef }
        open $fh, $file or return undef;
        $close = 1;
    }
    else
    {
        $fh = \*STDIN;
        $close = 0;
    }

    my $out = '';
    my $inc = 1048576;
    my $end =       0;
    my $read;
    while ( $read = read( $fh, $out, $inc, $end ) ) { $end += $read }
    close( $fh ) if $close;

    $out;
}


#-----------------------------------------------------------------------------
#  Previous, 50% slower fasta reader:
#-----------------------------------------------------------------------------
sub read_fasta_0
{
    my ( $fh, $name, $close, $unused ) = input_filehandle( $_[0] );
    $unused && die "Bad reference type (" . ref( $unused ) . ") passed to read_fasta\n";

    my @seqs = ();
    my ($id, $desc, $seq) = ("", "", "");

    while ( <$fh> ) {
        chomp;
        if (/^>\s*(\S+)(\s+(.*))?$/) {        #  new id
            if ($id && $seq) { push @seqs, [ $id, $desc, $seq ] }
            ($id, $desc, $seq) = ($1, $3 ? $3 : "", "");
        }
        else {
            tr/     0-9//d;
            $seq .= $_ ;
        }
    }
    close( $fh ) if $close;

    if ( $id && $seq ) { push @seqs, [ $id, $desc, $seq ] }
    return wantarray ? @seqs : \@seqs;
}


#-----------------------------------------------------------------------------
#  Read one fasta sequence at a time from a file.  This is half as fast a
#  read_fasta(), but can handle an arbitrarily large file.  State information
#  is retained in hashes, so any number of streams can be interlaced.
#
#      @entry = read_next_fasta_seq( \*FILEHANDLE )
#     \@entry = read_next_fasta_seq( \*FILEHANDLE )
#      @entry = read_next_fasta_seq(  $filename )
#     \@entry = read_next_fasta_seq(  $filename )
#      @entry = read_next_fasta_seq()                # \*STDIN
#     \@entry = read_next_fasta_seq()                # \*STDIN
#
#      @entry = ( $id, $description, $seq )
#
#  When reading at the end of file, () is returned.
#  With a filename, reading past this will reopen the file at the beginning.
#-----------------------------------------------------------------------------
#  Reading always overshoots, so save next id and description

{   #  Use bare block to scope the header hash

    my %next_header;
    my %file_handle;
    my %close_file;

    sub read_next_fasta_seq
    {
        $_[0] ||= \*STDIN;               #  Undefined $_[0] fails with use warn
        my $fh = $file_handle{ $_[0] };
        if ( ! $fh )
        {
            if ( ref $_[0] )
            {
                return () if ref $_[0] ne 'GLOB';
                $fh = $_[0];
            }
            else
            {
                my $file = $_[0];
                if    ( -f $file                       ) { $file = "<$file" }
                elsif (    $file =~ /^<(.*)$/ && -f $1 ) { }  # Explicit read
              # elsif (    $file =~ /\S\s*\|$/         ) { }  # Read from a pipe
                else                                     { return () }
                open $fh, $file or return ();
                $close_file{ $fh } = 1;
            }
            $file_handle{ $_[0] } = $fh;
        }

        my ( $id, $desc, $seq ) = ( undef, '', '' );
        if ( defined( $next_header{$fh} ) )
        {
            ( $id, $desc ) = parse_fasta_title( $next_header{$fh} );
        }
        else
        {
            $next_header{$fh} = '';
        }

        while ( <$fh> )
        {
            chomp;
            if ( /^>/ )        #  new id
            {
                $next_header{$fh} = $_;
                if ( defined($id) && $seq )
                {
                    return wantarray ? ($id, $desc, $seq) : [$id, $desc, $seq]
                }
                ( $id, $desc ) = parse_fasta_title( $next_header{$fh} );
                $seq = '';
            }
            else
            {
                tr/ \t\r//d;
                $seq .= $_;
            }
        }

        #  Done with file; there is no next header:

        delete $next_header{ $fh };

        #  Return last set of data:

        if ( defined($id) && $seq )
        {
            return wantarray ? ($id,$desc,$seq) : [$id,$desc,$seq]
        }

        #  Or close everything out (returning the empty list tells caller
        #  that we are done)

        if ( $close_file{ $fh } ) { close $fh; delete $close_file{ $fh } }
        delete $file_handle{ $_[0] };

        return ();
    }
}


#-----------------------------------------------------------------------------
#  Read a clustal alignment from a file.
#  Save the contents in a list of refs to arrays:  (id, description, seq)
#
#     @seq_entries = read_clustal_file( $filename )
#-----------------------------------------------------------------------------
sub read_clustal_file { read_clustal( @_ ) }


#-----------------------------------------------------------------------------
#  Read a clustal alignment.
#  Save the contents in a list of refs to arrays:  (id, description, seq)
#
#     @seq_entries = read_clustal( )              # STDIN
#     @seq_entries = read_clustal( \*FILEHANDLE )
#     @seq_entries = read_clustal(  $filename )
#-----------------------------------------------------------------------------
sub read_clustal {
    my ( $fh, undef, $close, $unused ) = input_filehandle( shift );
    $unused && die "Bad reference type (" . ref( $unused ) . ") passed to read_clustal_file\n";

    my ( %seq, @ids, $line );
    while ( defined( $line = <$fh> ) )
    {
        ( $line =~ /^[A-Za-z0-9]/ ) or next;
        chomp $line;
        my @flds = split /\s+/, $line;
        if ( @flds == 2 )
        {
            $seq{ $flds[0] } or push @ids, $flds[0];
            push @{ $seq{ $flds[0] } }, $flds[1];
        }
    }
    close( $fh ) if $close;

    map { [ $_, "", join( "", @{$seq{$_}} ) ] } @ids;
}


#-----------------------------------------------------------------------------
#  Parse a fasta file header to id and definition parts
#
#     ($id, $def) = parse_fasta_title( $title )
#     ($id, $def) = split_fasta_title( $title )
#-----------------------------------------------------------------------------
sub parse_fasta_title
{
    my $title = shift;
    chomp $title;

    return $title =~ /^>?\s*(\S+)(\s+(.*\S)?\s*)?$/ ? ( $1, $3 || '' )
         : $title =~ /^>/                           ? ( '', '' )
         :                                            ( undef, undef )
}

sub split_fasta_title { parse_fasta_title( @_ ) }


#-----------------------------------------------------------------------------
#  Helper function for defining an output filehandle:
#     filehandle is passed through
#     string is taken as file name to be openend
#     undef or "" defaults to STDOUT
#
#    ( \*FH, $close [, $file] ) = output_filehandle( $file );
#
#-----------------------------------------------------------------------------
sub output_filehandle
{
    my $file = shift;

    #  Null string or undef

    return ( \*STDOUT, 0 ) if ( ! defined( $file ) || ( $file eq "" ) );

    #  FILEHANDLE

    return ( $file, 0 ) if ( ref( $file ) eq "GLOB" );

    #  Some other kind of reference; return the unused value

    return ( \*STDOUT, 0, $file ) if ref( $file );

    #  File name

    my $fh;
    open( $fh, ">$file" ) || die "Could not open output $file\n";
    return ( $fh, 1 );
}


#-----------------------------------------------------------------------------
#  Legacy function for printing fasta sequence set:
#
#     print_seq_list_as_fasta( \*FILEHANDLE, @seq_entry_list );
#-----------------------------------------------------------------------------
sub print_seq_list_as_fasta { print_alignment_as_fasta( @_ ) }


#-----------------------------------------------------------------------------
#  Print list of sequence entries in fasta format.
#  Missing, undef or "" filename defaults to STDOUT.
#
#     print_alignment_as_fasta(                @seq_entry_list ); # STDOUT
#     print_alignment_as_fasta(               \@seq_entry_list ); # STDOUT
#     print_alignment_as_fasta( \*FILEHANDLE,  @seq_entry_list );
#     print_alignment_as_fasta( \*FILEHANDLE, \@seq_entry_list );
#     print_alignment_as_fasta(  $filename,    @seq_entry_list );
#     print_alignment_as_fasta(  $filename,   \@seq_entry_list );
#-----------------------------------------------------------------------------
sub print_alignment_as_fasta {
    my ( $fh, $close, $unused ) = output_filehandle( shift );
    ( unshift @_, $unused ) if $unused;

    ( ref( $_[0] ) eq "ARRAY" ) or confess "Bad sequence entry passed to print_alignment_as_fasta\n";

    #  Expand the sequence entry list if necessary:

    if ( ref( $_[0]->[0] ) eq "ARRAY" ) { @_ = @{ $_[0] } }

    foreach my $seq_ptr ( @_ ) { print_seq_as_fasta( $fh, @$seq_ptr ) }

    close( $fh ) if $close;
}


#-----------------------------------------------------------------------------
#  Print list of sequence entries in phylip format.
#  Missing, undef or "" filename defaults to STDOUT.
#
#     print_alignment_as_phylip(                @seq_entry_list ); # STDOUT
#     print_alignment_as_phylip(               \@seq_entry_list ); # STDOUT
#     print_alignment_as_phylip( \*FILEHANDLE,  @seq_entry_list );
#     print_alignment_as_phylip( \*FILEHANDLE, \@seq_entry_list );
#     print_alignment_as_phylip(  $filename,    @seq_entry_list );
#     print_alignment_as_phylip(  $filename,   \@seq_entry_list );
#-----------------------------------------------------------------------------
sub print_alignment_as_phylip {
    my ( $fh, $close, $unused ) = output_filehandle( shift );
    ( unshift @_, $unused ) if $unused;

    ( ref( $_[0] ) eq "ARRAY" ) or die die "Bad sequence entry passed to print_alignment_as_phylip\n";

    my @seq_list = ( ref( $_[0]->[0] ) eq "ARRAY" ) ? @{ $_[0] } : @_;

    my ( %id2, %used );
    my $maxlen = 0;
    foreach ( @seq_list )
    {
        my ( $id, undef, $seq ) = @$_;

        #  Need a name that is unique within 10 characters

        my $id2 = substr( $id, 0, 10 );
        $id2 =~ s/_/ /g;  # PHYLIP sequence files accept spaces
        my $n = "0";
        while ( $used{ $id2 } )
        {
            $n++;
            $id2 = substr( $id, 0, 10 - length( $n ) ) . $n;
        }
        $used{ $id2 } = 1;
        $id2{ $id } = $id2;

                #  Prepare to pad sequences (should not be necessary, but ...)

        my $len = length( $seq );
        $maxlen = $len if ( $len > $maxlen );
    }

    my $nseq = @seq_list;
    print $fh "$nseq  $maxlen\n";
    foreach ( @seq_list )
    {
        my ( $id, undef, $seq ) = @$_;
        my $len = length( $seq );
        printf $fh "%-10s  %s%s\n", $id2{ $id },
                                    $seq,
                                    $len<$maxlen ? ("?" x ($maxlen-$len)) : "";
    }

    close( $fh ) if $close;
}


#-----------------------------------------------------------------------------
#  Print list of sequence entries in nexus format.
#  Missing, undef or "" filename defaults to STDOUT.
#
#     print_alignment_as_nexus(               [ \%label_hash, ]  @seq_entry_list );
#     print_alignment_as_nexus(               [ \%label_hash, ] \@seq_entry_list );
#     print_alignment_as_nexus( \*FILEHANDLE, [ \%label_hash, ]  @seq_entry_list );
#     print_alignment_as_nexus( \*FILEHANDLE, [ \%label_hash, ] \@seq_entry_list );
#     print_alignment_as_nexus(  $filename,   [ \%label_hash, ]  @seq_entry_list );
#     print_alignment_as_nexus(  $filename,   [ \%label_hash, ] \@seq_entry_list );
#-----------------------------------------------------------------------------
sub print_alignment_as_nexus {
    my ( $fh, $close, $unused ) = output_filehandle( shift );
    ( unshift @_, $unused ) if $unused;

    my $lbls = ( ref( $_[0] ) eq "HASH" ) ? shift : undef;

    ( ref( $_[0] ) eq "ARRAY" ) or die "Bad sequence entry passed to print_alignment_as_nexus\n";

    my @seq_list = ( ref( $_[0]->[0] ) eq "ARRAY" ) ? @{ $_[0] } : @_;

    my %id2;
    my ( $maxidlen, $maxseqlen ) = ( 0, 0 );
    my ( $n1, $n2, $nt, $nu ) = ( 0, 0, 0, 0 );
    foreach ( @seq_list )
    {
        my ( $id, undef, $seq ) = @$_;
        my $id2 = $lbls ? ( $lbls->{ $id } || $id ) : $id;
        if ( $id2 !~ /^[-+.0-9A-Za-z~_|]+$/ )
        {
                $id2 =~ s/'/''/g;
                $id2 = qq('$id2');
            }
        $id2{ $id } = $id2;
        my $idlen = length( $id2 );
        $maxidlen = $idlen if ( $idlen > $maxidlen );

        my $seqlen = length( $seq );
        $maxseqlen = $seqlen if ( $seqlen > $maxseqlen );

        $nt += $seq =~ tr/Tt//d;
        $nu += $seq =~ tr/Uu//d;
        $n1 += $seq =~ tr/ACGNacgn//d;
        $n2 += $seq =~ tr/A-Za-z//d;
    }

    my $nseq = @seq_list;
    my $type = ( $n1 < 2 * $n2 ) ?  'protein' : ($nt>$nu) ? 'DNA' : 'RNA';

    print $fh <<"END_HEAD";
#NEXUS

BEGIN Data;
    Dimensions
        NTax=$nseq
        NChar=$maxseqlen
        ;
    Format
        DataType=$type
        Gap=-
        Missing=?
        ;
    Matrix

END_HEAD

    foreach ( @seq_list )
    {
        my ( $id, undef, $seq ) = @$_;
        my $len = length( $seq );
        printf  $fh  "%-${maxidlen}s  %s%s\n",
                     $id2{ $id },
                     $seq,
                     $len<$maxseqlen ? ("?" x ($maxseqlen-$len)) : "";
    }

    print $fh <<"END_TAIL";
;
END;
END_TAIL

    close( $fh ) if $close;
}


#-----------------------------------------------------------------------------
#  Print one sequence in fasta format to an open file.
#
#     print_seq_as_fasta( \*FILEHANDLE, $id, $desc, $seq );
#     print_seq_as_fasta(               $id, $desc, $seq );
#     print_seq_as_fasta( \*FILEHANDLE, $id,        $seq );
#     print_seq_as_fasta(               $id,        $seq );
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  print_seq_as_fasta() is meant more as a internal support routine than an
#  external interface.  To print a single sequence to a named file use:
#
#     print_alignment_as_fasta( $filename, [ $id, $desc, $seq ] );
#     print_alignment_as_fasta( $filename, [ $id,        $seq ] );
#-----------------------------------------------------------------------------
sub print_seq_as_fasta
{
    my $fh = ( ref $_[0] eq 'GLOB' ) ? shift : \*STDOUT;
    return if ( @_ < 2 ) || ( @_ > 3 ) || ! ( defined $_[0] && defined $_[-1] );
    #  Print header line
    print $fh  ( @_ == 3 && defined $_[1] && $_[1] =~ /\S/ ) ? ">$_[0] $_[1]\n" : ">$_[0]\n";
    #  Print sequence, 60 chars per line
    print $fh  join( "\n", $_[-1] =~ m/.{1,60}/g ), "\n";
}


#-----------------------------------------------------------------------------
#  Print one sequence in GenBank flat file format:
#
#     print_gb_locus( \*FILEHANDLE, $locus, $def, $accession, $seq )
#-----------------------------------------------------------------------------
sub print_gb_locus {
    my ($fh, $loc, $def, $acc, $seq) = @_;
    my ($len, $i, $imax);
    my $istep = 10;

    $len = length($seq);
    printf $fh  "LOCUS       %-10s%7d bp\n", substr($loc,0,10), $len;
    print  $fh  "DEFINITION  " . substr(wrap_text($def,80,12), 12) . "\n";
    if ($acc) { print  $fh  "ACCESSION   $acc\n" }
    print  $fh "ORIGIN\n";

    for ($i = 1; $i <= $len; ) {
        printf $fh "%9d", $i;
        $imax = $i + 59; if ($imax > $len) { $imax = $len }
        for ( ; $i <= $imax; $i += $istep) {
            print $fh " " . substr($seq, $i-1, $istep);
        }
        print $fh "\n";
        }
    print $fh "//\n";
}


#-----------------------------------------------------------------------------
#  Return a string with text wrapped to defined line lengths:
#
#     $wrapped_text = wrap_text( $str )                  # default len   =  80
#     $wrapped_text = wrap_text( $str, $len )            # default ind   =   0
#     $wrapped_text = wrap_text( $str, $len, $indent )   # default ind_n = ind
#     $wrapped_text = wrap_text( $str, $len, $indent_1, $indent_n )
#-----------------------------------------------------------------------------
sub wrap_text {
    my ($str, $len, $ind, $indn) = @_;

    defined($str)  || die "wrap_text called without a string\n";
    defined($len)  || ($len  =   80);
    defined($ind)  || ($ind  =    0);
    ($ind  < $len) || die "wrap error: indent greater than line length\n";
    defined($indn) || ($indn = $ind);
    ($indn < $len) || die "wrap error: indent_n greater than line length\n";

    $str =~ s/\s+$//;
    $str =~ s/^\s+//;
    my ($maxchr, $maxchr1);
    my (@lines) = ();

    while ($str) {
        $maxchr1 = ($maxchr = $len - $ind) - 1;
        if ($maxchr >= length($str)) {
            push @lines, (" " x $ind) . $str;
            last;
        }
        elsif ($str =~ /^(.{0,$maxchr1}\S)\s+(\S.*)$/) { # no expr in {}
            push @lines, (" " x $ind) . $1;
            $str = $2;
        }
        elsif ($str =~ /^(.{0,$maxchr1}-)(.*)$/) {
            push @lines, (" " x $ind) . $1;
            $str = $2;
        }
        else {
            push @lines, (" " x $ind) . substr($str, 0, $maxchr);
            $str = substr($str, $maxchr);
        }
        $ind = $indn;
    }

    return join("\n", @lines);
}


#-----------------------------------------------------------------------------
#  Build an index from seq_id to pointer to sequence entry: (id, desc, seq)
#
#     my \%seq_ind  = index_seq_list(  @seq_list );
#     my \%seq_ind  = index_seq_list( \@seq_list );
#
#  Usage example:
#
#  my  @seq_list   = read_fasta_seqs(\*STDIN);  # list of pointers to entries
#  my \%seq_ind    = index_seq_list(@seq_list); # hash from names to pointers
#  my  @chosen_seq = @{%seq_ind{"contig1"}};    # extract one entry
#
#-----------------------------------------------------------------------------
sub index_seq_list {
    ( ref( $_[0] )      ne 'ARRAY' ) ? {}
  : ( ref( $_[0]->[0] ) ne 'ARRAY' ) ? { map { $_->[0] => $_ } @_ }
  :                                    { map { $_->[0] => $_ } @{ $_[0] } }
}


#-----------------------------------------------------------------------------
#  Three routines to access all or part of sequence entry by id:
#
#     @seq_entry = seq_entry_by_id( \%seq_index, $seq_id );
#    \@seq_entry = seq_entry_by_id( \%seq_index, $seq_id );
#     $seq_desc  = seq_desc_by_id(  \%seq_index, $seq_id );
#     $seq       = seq_data_by_id(  \%seq_index, $seq_id );
#
#-----------------------------------------------------------------------------
sub seq_entry_by_id {
    (my $ind_ref = shift)  || die "No index supplied to seq_entry_by_id\n";
    (my $id      = shift)  || die "No id supplied to seq_entry_by_id\n";
    return wantarray ? @{ $ind_ref->{$id} } : $ind_ref->{$id};
}


sub seq_desc_by_id {
    (my $ind_ref = shift)  || die "No index supplied to seq_desc_by_id\n";
    (my $id      = shift)  || die "No id supplied to seq_desc_by_id\n";
    return ${ $ind_ref->{$id} }[1];
}


sub seq_data_by_id {
    (my $ind_ref = shift)  || die "No index supplied to seq_data_by_id\n";
    (my $id      = shift)  || die "No id supplied to seq_data_by_id\n";
    return ${ $ind_ref->{$id} }[2];
}


#-----------------------------------------------------------------------------
#  Remove columns of alignment gaps from sequences:
#
#   @packed_seqs = pack_alignment(  @seqs )
#   @packed_seqs = pack_alignment( \@seqs )
#  \@packed_seqs = pack_alignment(  @seqs )
#  \@packed_seqs = pack_alignment( \@seqs )
#
#  Gap characters are defined below as '-', '~', '.', and ' '.
#-----------------------------------------------------------------------------
sub pack_alignment
{
    $_[0] && ( ref( $_[0] ) eq 'ARRAY' ) && @{$_[0]} && defined( $_[0]->[0] )
        or return ();

    my @seqs = ( ref( $_[0]->[0] ) eq 'ARRAY' ) ? @{$_[0] } : @_;
    @seqs or return wantarray ? () : [];

    my $mask  = gap_mask( $seqs[0]->[2] );
    foreach ( @seqs[ 1 .. (@seqs-1) ] )
    {
        $mask |= gap_mask( $_->[2] );
    }

    my $seq;
    my @seqs2 = map { $seq = $_->[2] & $mask;
                      $seq =~ tr/\000//d;
                      [ $_->[0], $_->[1], $seq ]
                    }
                @seqs;

    wantarray ? @seqs2 : \@seqs2;
}


#-------------------------------------------------------------------------------
#  Produce a packing mask for columns of alignment gaps in sequences.  Gap
#  columns are 0x00 characters, and all others are 0xFF.
#
#   $mask = alignment_gap_mask(  @seqs )
#   $mask = alignment_gap_mask( \@seqs )
#
#-------------------------------------------------------------------------------
sub alignment_gap_mask
{
    $_[0] && ( ref( $_[0] ) eq 'ARRAY' ) && @{$_[0]} && defined( $_[0]->[0] )
        or return undef;

    my @seqs = ( ref( $_[0]->[0] ) eq 'ARRAY' ) ? @{$_[0] } : @_;
    @seqs or return undef;

    my $mask = gap_mask( $seqs[0]->[2] );
    foreach ( @seqs[ 1 .. (@seqs-1) ] ) { $mask |= gap_mask( $_->[2] ) }

    $mask;
}


#-------------------------------------------------------------------------------
#  Pack a sequence alignment according to a mask, removing positions where
#  mask has 0x00 (or '-') characters
#
#   @packed = pack_alignment_by_mask( $mask,  @align )
#   @packed = pack_alignment_by_mask( $mask, \@align )
#  \@packed = pack_alignment_by_mask( $mask,  @align )
#  \@packed = pack_alignment_by_mask( $mask, \@align )
#
#-------------------------------------------------------------------------------
sub pack_alignment_by_mask
{
    my $mask = shift;
    defined $mask && ! ref( $mask ) && length( $mask )
        or return ();
    $mask =~ tr/-/\000/;      # Allow '-' as a column to be removed
    $mask =~ tr/\000/\377/c;  # Make sure that everything not null is 0xFF
 
    $_[0] && ( ref( $_[0] ) eq 'ARRAY' ) && @{$_[0]} && defined( $_[0]->[0] )
        or return ();
    my @seqs = ( ref( $_[0]->[0] ) eq 'ARRAY' ) ? @{$_[0] } : @_;

    my $seq;
    my @seqs2 = map { $seq = $_->[2] & $mask;     # Apply mask to sequence
                      $seq =~ tr/\000//d;         # Delete null characters
                      [ $_->[0], $_->[1], $seq ]  # Rebuild sequence entries
                    }
                @seqs;

    wantarray ? @seqs2 : \@seqs2;
}


#-------------------------------------------------------------------------------
#  Weight a sequence alignment according to a mask of digits, 0-9.
#
#   @packed = weight_alignment_by_mask( $mask,  @align )
#   @packed = weight_alignment_by_mask( $mask, \@align )
#  \@packed = weight_alignment_by_mask( $mask,  @align )
#  \@packed = weight_alignment_by_mask( $mask, \@align )
#
#-------------------------------------------------------------------------------
sub weight_alignment_by_mask
{
    my $mask = shift;
    defined $mask && ! ref( $mask ) && length( $mask )
        or return ();
 
    $_[0] && ( ref( $_[0] ) eq 'ARRAY' ) && @{$_[0]} && defined( $_[0]->[0] )
        or return ();
    my @seqs = ( ref( $_[0]->[0] ) eq 'ARRAY' ) ? @{$_[0] } : @_;

    my @seqs2 = map { [ $_->[0], $_->[1], weight_seq_by_mask_0( $mask, $_->[2] ) ] } @seqs;

    wantarray ? @seqs2 : \@seqs2;
}


#
#  Assume data are valid
#
sub weight_seq_by_mask_0
{
    my ( $mask, $seq ) = @_;

    #  Remove 0 weight columns, which is fast and easy:
    my $m0 = $mask;
    $m0 =~ tr/123456789/\377/;
    $m0 =~ tr/\377/\000/c;
    ( $mask &= $m0 ) =~ tr/\000//d;
    ( $seq  &= $m0 ) =~ tr/\000//d;

    #  If all remaining cols are weight 1, we are done:
    return $seq if $mask =~ /^1*$/;

    my @seq;
    for ( my $i = 0; $i < length( $mask ); $i++ )
    {
        push @seq, substr( $seq, $i, 1 ) x substr( $mask, $i, 1 );
    }

    join( '', @seq );
}


#-----------------------------------------------------------------------------
#  Make a mask in which gap characters ('-', '~', '.', and ' ') are converted
#  to 0x00, and all other characters to 0xFF.
#
#      $mask = gap_mask( $seq )
#
#-----------------------------------------------------------------------------
sub gap_mask
{
    my $mask = shift;
    defined $mask or return '';

    $mask =~ tr/-~. /\000/;    #  Gap characters (might be extended)
    $mask =~ tr/\000/\377/c;   #  Non-gap characters
    $mask;
}


#===============================================================================
#  Expand sequences with the local gap character in a manner that reverses the
#  pack by mask function.
#
#      $expanded = expand_sequence_by_mask( $seq, $mask )
#
#  The columns to be added can be marked by '-' or "\000" in the mask. 
#
#  Code note:
#
#  The function attempts to match the local gap character in the sequence.
#  $c0 and $c1 are the previous and next characters in the sequence being
#  expanded.  (($c0,$c1)=($c1,shift @s1))[0] updates the values and evaluates
#  to what was the next character, and becomes the new previous character.
#  The following really does print w, x, y, and z, one per line:
#
#     ( $a, $b, @c ) = ("", split //, "wxyz");
#     while ( defined $b ) { print( (($a,$b)=($b,shift @c))[0], "\n" ) }
#===============================================================================

sub expand_sequence_by_mask
{
    my ( $seq, $mask ) = @_;

    $mask =~ tr/-/\000/;  #  Allow hyphen or null in mask at added positions.
    my ( $c0, $c1, @s1 ) = ( '', split( //, $seq ), '' );

    join '', map { $_  ne "\000"            ? (($c0,$c1)=($c1,shift @s1))[0]
                 : $c0 eq '~' || $c1 eq '~' ? '~'
                 : $c0 eq '.' || $c1 eq '.' ? '.'
                 :                            '-'
                 }
             split //, $mask;
}


#-----------------------------------------------------------------------------
#  Remove all alignment gaps from sequences:
#
#   @packed_seqs = pack_sequences(  @seqs )  # Also handles single sequence
#   @packed_seqs = pack_sequences( \@seqs )
#  \@packed_seqs = pack_sequences(  @seqs )
#  \@packed_seqs = pack_sequences( \@seqs )
#
#-----------------------------------------------------------------------------
sub pack_sequences
{
    $_[0] && ( ref( $_[0] ) eq 'ARRAY' ) && @{$_[0]} && defined( $_[0]->[0] )
        or return ();

    my @seqs = ( ref( $_[0]->[0] ) eq 'ARRAY' ) ? @{$_[0] } : @_;

    my @seqs2 = map { [ $_->[0], $_->[1], pack_seq( $_->[2] ) ] } @seqs;

    wantarray ? @seqs2 : \@seqs2;
}


#-----------------------------------------------------------------------------
#  Some simple sequence manipulations:
#
#     @entry  = subseq_DNA_entry( @seq_entry, $from, $to [, $fix_id] );
#     @entry  = subseq_RNA_entry( @seq_entry, $from, $to [, $fix_id] );
#     @entry  = complement_DNA_entry( @seq_entry [, $fix_id] );
#     @entry  = complement_RNA_entry( @seq_entry [, $fix_id] );
#     $DNAseq = complement_DNA_seq( $NA_seq );
#     $RNAseq = complement_RNA_seq( $NA_seq );
#     $DNAseq = to_DNA_seq( $NA_seq );
#     $RNAseq = to_RNA_seq( $NA_seq );
#
#-----------------------------------------------------------------------------

sub subseq_DNA_entry {
    my ($id, $desc, @rest) = @_;
    wantarray || die "subseq_DNA_entry requires array context\n";

    my $seq;
    ($id, $seq) = subseq_nt(1, $id, @rest);  # 1 is for DNA, not RNA
    return ($id, $desc, $seq);
}


sub subseq_RNA_entry {
    my ($id, $desc, @rest) = @_;
    wantarray || die "subseq_RNA_entry requires array context\n";

    my $seq;
    ($id, $seq) = subseq_nt(0, $id, @rest);  # 0 is for not DNA, i.e., RNA
    return ($id, $desc, $seq);
}


sub subseq_nt {
    my ($DNA, $id, $seq, $from, $to, $fix_id) = @_;
    $fix_id ||= 0;     #  fix undef value

    my $len   = length($seq);
    if (         ( $from eq '$' ) || ( $from eq "" ) ) { $from = $len }
    if (! $to || ( $to   eq '$' ) || ( $to   eq "" ) ) { $to   = $len }

    my $left  = ( $from < $to ) ? $from : $to;
    my $right = ( $from < $to ) ? $to   : $from;
    if ( ( $right < 1 ) || ( $left > $len ) ) { return ($id, "") }
    if ( $right > $len ) { $right = $len }
    if ( $left  < 1    ) { $left  = 1 }

    $seq = substr($seq, $left-1, $right-$left+1);
    if ( $from > $to ) {
        $seq = reverse $seq;
        if ( $DNA ) {
            $seq =~ tr[ACGTUKMRSWYBDHVNacgtukmrswybdhvn]
                      [TGCAAMKYSWRVHDBNtgcaamkyswrvhdbn];
        }
        else {
            $seq =~ tr[ACGTUKMRSWYBDHVNacgtukmrswybdhvn]
                      [UGCAAMKYSWRVHDBNugcaamkyswrvhdbn];
        }
    }

    if ( $fix_id ) {
        if ( ( $id =~ s/_(\d+)_(\d+)$// )
          && ( abs($2-$1)+1 == $len ) ) {
            if ( $1 <= $2 ) { $from += $1 - 1;         $to += $1 - 1 }
            else            { $from  = $1 + 1 - $from; $to  = $1 + 1 - $to }
        }
        $id .= "_" . $from . "_" . $to;
    }

    return ($id, $seq);
}


sub DNA_subseq
{
    my ( $seq, $from, $to ) = @_;

    my $len = ref( $seq ) eq 'SCALAR' ? length( $$seq )
                                      : length(  $seq );
    if ( ( $from eq '$' ) || ( $from eq "" ) ) { $from = $len }
    if ( ( $to   eq '$' ) || ( ! $to       ) ) { $to   = $len }

    my $left  = ( $from < $to ) ? $from : $to;
    my $right = ( $from < $to ) ? $to   : $from;
    if ( ( $right < 1 ) || ( $left > $len ) ) { return "" }
    if ( $right > $len ) { $right = $len }
    if ( $left  < 1    ) { $left  =    1 }

    my $subseq = ref( $seq ) eq 'SCALAR' ? substr( $$seq, $left-1, $right-$left+1 )
                                         : substr(  $seq, $left-1, $right-$left+1 );

    if ( $from > $to )
    {
        $subseq = reverse $subseq;
        $subseq =~ tr[ACGTUKMRSWYBDHVNacgtukmrswybdhvn]
                     [TGCAAMKYSWRVHDBNtgcaamkyswrvhdbn];
    }

    $subseq
}


sub RNA_subseq
{
    my ( $seq, $from, $to ) = @_;

    my $len = ref( $seq ) eq 'SCALAR' ? length( $$seq )
                                      : length(  $seq );
    if ( ( $from eq '$' ) || ( $from eq "" ) ) { $from = $len }
    if ( ( $to   eq '$' ) || ( ! $to       ) ) { $to   = $len }

    my $left  = ( $from < $to ) ? $from : $to;
    my $right = ( $from < $to ) ? $to   : $from;
    if ( ( $right < 1 ) || ( $left > $len ) ) { return "" }
    if ( $right > $len ) { $right = $len }
    if ( $left  < 1    ) { $left  =    1 }

    my $subseq = ref( $seq ) eq 'SCALAR' ? substr( $$seq, $left-1, $right-$left+1 )
                                         : substr(  $seq, $left-1, $right-$left+1 );

    if ( $from > $to )
    {
        $subseq = reverse $subseq;
        $subseq =~ tr[ACGTUKMRSWYBDHVNacgtukmrswybdhvn]
                     [UGCAAMKYSWRVHDBNugcaamkyswrvhdbn];
    }

    $subseq
}


sub complement_DNA_entry {
    my ($id, $desc, $seq, $fix_id) = @_;
    $fix_id ||= 0;     #  fix undef values

    wantarray || die "complement_DNA_entry requires array context\n";
    $seq = reverse $seq;
    $seq =~ tr[ACGTUKMRSWYBDHVNacgtukmrswybdhvn]
              [TGCAAMKYSWRVHDBNtgcaamkyswrvhdbn];
    if ($fix_id) {
        if ($id =~ s/_(\d+)_(\d+)$//) {
            $id .= "_" . $2 . "_" . $1;
        }
        else {
            $id .= "_" . length($seq) . "_1";
        }
    }

    return ($id, $desc, $seq);
}


sub complement_RNA_entry {
    my ($id, $desc, $seq, $fix_id) = @_;
    $fix_id ||= 0;     #  fix undef values

    wantarray || die "complement_DNA_entry requires array context\n";
    $seq = reverse $seq;
    $seq =~ tr[ACGTUKMRSWYBDHVNacgtukmrswybdhvn]
              [UGCAAMKYSWRVHDBNugcaamkyswrvhdbn];
    if ($fix_id) {
        if ($id =~ s/_(\d+)_(\d+)$//) {
            $id .= "_" . $2 . "_" . $1;
        }
        else {
            $id .= "_" . length($seq) . "_1";
        }
    }

    return ($id, $desc, $seq);
}


sub complement_DNA_seq {
    my $seq = reverse shift;
    $seq =~ tr[ACGTUKMRSWYBDHVNacgtukmrswybdhvn]
              [TGCAAMKYSWRVHDBNtgcaamkyswrvhdbn];
    return $seq;
}


sub complement_RNA_seq {
    my $seq = reverse shift;
    $seq =~ tr[ACGTUKMRSWYBDHVNacgtukmrswybdhvn]
              [UGCAAMKYSWRVHDBNugcaamkyswrvhdbn];
    return $seq;
}


sub to_DNA_seq {
    my $seq = shift;
    $seq =~ tr/Uu/Tt/;
    return $seq;
}


sub to_RNA_seq {
    my $seq = shift;
    $seq =~ tr/Tt/Uu/;
    return $seq;
}


sub pack_seq {
    my $seq = shift;
    $seq =~ tr/A-Za-z*//cd;
    $seq;
}


sub clean_ae_sequence {
    local $_ = shift;
    $_ = to7bit($_);
    s/\+/1/g;
    s/[^0-9A-IK-NP-Za-ik-np-z~.-]/-/g;
    return $_;
}


sub to7bit {
    local $_ = shift;
    my ($o, $c);
    while (/\\([0-3][0-7][0-7])/) {
        $o = oct($1) % 128;
        $c = sprintf("%c", $o);
        s/\\$1/$c/g;
    }
    return $_;
}


sub to8bit {
    local $_ = shift;
    my ($o, $c);
    while (/\\([0-3][0-7][0-7])/) {
        $o = oct($1);
        $c = sprintf("%c", $o);
        s/\\$1/$c/g;
    }
    return $_;
}



#-----------------------------------------------------------------------------
#  Translate nucleotides to one letter protein:
#
#     $seq = translate_seq( $seq [, $met_start] )
#     $aa  = translate_codon( $triplet )
#     $aa  = translate_DNA_codon( $triplet )     # Does not rely on DNA
#     $aa  = translate_uc_DNA_codon( $triplet )  # Does not rely on uc or DNA
#
#  User-supplied genetic code must be upper case index and match the
#  DNA versus RNA type of sequence
#
#     $seq = translate_seq_with_user_code( $seq, $gen_code_hash [, $met_start] )
#
#-----------------------------------------------------------------------------

@aa_1_letter_order = qw( A C D E F G H I K L M N P Q R S T V W Y );  # Alpha by 1 letter
@aa_3_letter_order = qw( A R N D C Q E G H I L K M F P S T W Y V );  # PAM matrix order
@aa_n_codon_order  = qw( L R S A G P T V I C D E F H K N Q Y M W );  

%genetic_code = (

    # DNA version

    TTT => 'F',  TCT => 'S',  TAT => 'Y',  TGT => 'C',
    TTC => 'F',  TCC => 'S',  TAC => 'Y',  TGC => 'C',
    TTA => 'L',  TCA => 'S',  TAA => '*',  TGA => '*',
    TTG => 'L',  TCG => 'S',  TAG => '*',  TGG => 'W',
    CTT => 'L',  CCT => 'P',  CAT => 'H',  CGT => 'R',
    CTC => 'L',  CCC => 'P',  CAC => 'H',  CGC => 'R',
    CTA => 'L',  CCA => 'P',  CAA => 'Q',  CGA => 'R',
    CTG => 'L',  CCG => 'P',  CAG => 'Q',  CGG => 'R',
    ATT => 'I',  ACT => 'T',  AAT => 'N',  AGT => 'S',
    ATC => 'I',  ACC => 'T',  AAC => 'N',  AGC => 'S',
    ATA => 'I',  ACA => 'T',  AAA => 'K',  AGA => 'R',
    ATG => 'M',  ACG => 'T',  AAG => 'K',  AGG => 'R',
    GTT => 'V',  GCT => 'A',  GAT => 'D',  GGT => 'G',
    GTC => 'V',  GCC => 'A',  GAC => 'D',  GGC => 'G',
    GTA => 'V',  GCA => 'A',  GAA => 'E',  GGA => 'G',
    GTG => 'V',  GCG => 'A',  GAG => 'E',  GGG => 'G',

    #  The following ambiguous encodings are not necessary,  but
    #  speed up the processing of some ambiguous triplets:

    TTY => 'F',  TCY => 'S',  TAY => 'Y',  TGY => 'C',
    TTR => 'L',  TCR => 'S',  TAR => '*',
                 TCN => 'S',
    CTY => 'L',  CCY => 'P',  CAY => 'H',  CGY => 'R',
    CTR => 'L',  CCR => 'P',  CAR => 'Q',  CGR => 'R',
    CTN => 'L',  CCN => 'P',               CGN => 'R',
    ATY => 'I',  ACY => 'T',  AAY => 'N',  AGY => 'S',
                 ACR => 'T',  AAR => 'K',  AGR => 'R',
                 ACN => 'T',
    GTY => 'V',  GCY => 'A',  GAY => 'D',  GGY => 'G',
    GTR => 'V',  GCR => 'A',  GAR => 'E',  GGR => 'G',
    GTN => 'V',  GCN => 'A',               GGN => 'G'
);

#  Add RNA by construction:

foreach ( grep { /T/ } keys %genetic_code )
{
    my $codon = $_;
    $codon =~ s/T/U/g;
    $genetic_code{ $codon } = lc $genetic_code{ $_ }
}

#  Add lower case by construction:

foreach ( keys %genetic_code )
{
    $genetic_code{ lc $_ } = lc $genetic_code{ $_ }
}


#  Construct the genetic code with selenocysteine by difference:

%genetic_code_with_U = %genetic_code;
$genetic_code_with_U{ TGA } = 'U';
$genetic_code_with_U{ tga } = 'u';
$genetic_code_with_U{ UGA } = 'U';
$genetic_code_with_U{ uga } = 'u';


%amino_acid_codons_DNA = (
         L  => [ qw( TTA TTG CTA CTG CTT CTC ) ],
         R  => [ qw( AGA AGG CGA CGG CGT CGC ) ],
         S  => [ qw( AGT AGC TCA TCG TCT TCC ) ],
         A  => [ qw( GCA GCG GCT GCC ) ],
         G  => [ qw( GGA GGG GGT GGC ) ],
         P  => [ qw( CCA CCG CCT CCC ) ],
         T  => [ qw( ACA ACG ACT ACC ) ],
         V  => [ qw( GTA GTG GTT GTC ) ],
         I  => [ qw( ATA ATT ATC ) ],
         C  => [ qw( TGT TGC ) ],
         D  => [ qw( GAT GAC ) ],
         E  => [ qw( GAA GAG ) ],
         F  => [ qw( TTT TTC ) ],
         H  => [ qw( CAT CAC ) ],
         K  => [ qw( AAA AAG ) ],
         N  => [ qw( AAT AAC ) ],
         Q  => [ qw( CAA CAG ) ],
         Y  => [ qw( TAT TAC ) ],
         M  => [ qw( ATG ) ],
         U  => [ qw( TGA ) ],
         W  => [ qw( TGG ) ],

         l  => [ qw( tta ttg cta ctg ctt ctc ) ],
         r  => [ qw( aga agg cga cgg cgt cgc ) ],
         s  => [ qw( agt agc tca tcg tct tcc ) ],
         a  => [ qw( gca gcg gct gcc ) ],
         g  => [ qw( gga ggg ggt ggc ) ],
         p  => [ qw( cca ccg cct ccc ) ],
         t  => [ qw( aca acg act acc ) ],
         v  => [ qw( gta gtg gtt gtc ) ],
         i  => [ qw( ata att atc ) ],
         c  => [ qw( tgt tgc ) ],
         d  => [ qw( gat gac ) ],
         e  => [ qw( gaa gag ) ],
         f  => [ qw( ttt ttc ) ],
         h  => [ qw( cat cac ) ],
         k  => [ qw( aaa aag ) ],
         n  => [ qw( aat aac ) ],
         q  => [ qw( caa cag ) ],
         y  => [ qw( tat tac ) ],
         m  => [ qw( atg ) ],
         u  => [ qw( tga ) ],
         w  => [ qw( tgg ) ],

        '*' => [ qw( TAA TAG TGA ) ]
);



%amino_acid_codons_RNA = (
         L  => [ qw( UUA UUG CUA CUG CUU CUC ) ],
         R  => [ qw( AGA AGG CGA CGG CGU CGC ) ],
         S  => [ qw( AGU AGC UCA UCG UCU UCC ) ],
         A  => [ qw( GCA GCG GCU GCC ) ],
         G  => [ qw( GGA GGG GGU GGC ) ],
         P  => [ qw( CCA CCG CCU CCC ) ],
         T  => [ qw( ACA ACG ACU ACC ) ],
         V  => [ qw( GUA GUG GUU GUC ) ],
         B  => [ qw( GAU GAC AAU AAC ) ],
         Z  => [ qw( GAA GAG CAA CAG ) ],
         I  => [ qw( AUA AUU AUC ) ],
         C  => [ qw( UGU UGC ) ],
         D  => [ qw( GAU GAC ) ],
         E  => [ qw( GAA GAG ) ],
         F  => [ qw( UUU UUC ) ],
         H  => [ qw( CAU CAC ) ],
         K  => [ qw( AAA AAG ) ],
         N  => [ qw( AAU AAC ) ],
         Q  => [ qw( CAA CAG ) ],
         Y  => [ qw( UAU UAC ) ],
         M  => [ qw( AUG ) ],
         U  => [ qw( UGA ) ],
         W  => [ qw( UGG ) ],

         l  => [ qw( uua uug cua cug cuu cuc ) ],
         r  => [ qw( aga agg cga cgg cgu cgc ) ],
         s  => [ qw( agu agc uca ucg ucu ucc ) ],
         a  => [ qw( gca gcg gcu gcc ) ],
         g  => [ qw( gga ggg ggu ggc ) ],
         p  => [ qw( cca ccg ccu ccc ) ],
         t  => [ qw( aca acg acu acc ) ],
         v  => [ qw( gua gug guu guc ) ],
         b  => [ qw( gau gac aau aac ) ],
         z  => [ qw( gaa gag caa cag ) ],
         i  => [ qw( aua auu auc ) ],
         c  => [ qw( ugu ugc ) ],
         d  => [ qw( gau gac ) ],
         e  => [ qw( gaa gag ) ],
         f  => [ qw( uuu uuc ) ],
         h  => [ qw( cau cac ) ],
         k  => [ qw( aaa aag ) ],
         n  => [ qw( aau aac ) ],
         q  => [ qw( caa cag ) ],
         y  => [ qw( uau uac ) ],
         m  => [ qw( aug ) ],
         u  => [ qw( uga ) ],
         w  => [ qw( ugg ) ],

        '*' => [ qw( UAA UAG UGA ) ]
);


%n_codon_for_aa = map {
    $_ => scalar @{ $amino_acid_codons_DNA{ $_ } }
    } keys %amino_acid_codons_DNA;


%reverse_genetic_code_DNA = (
         A  => "GCN",  a  => "gcn",
         C  => "TGY",  c  => "tgy",
         D  => "GAY",  d  => "gay",
         E  => "GAR",  e  => "gar",
         F  => "TTY",  f  => "tty",
         G  => "GGN",  g  => "ggn",
         H  => "CAY",  h  => "cay",
         I  => "ATH",  i  => "ath",
         K  => "AAR",  k  => "aar",
         L  => "YTN",  l  => "ytn",
         M  => "ATG",  m  => "atg",
         N  => "AAY",  n  => "aay",
         P  => "CCN",  p  => "ccn",
         Q  => "CAR",  q  => "car",
         R  => "MGN",  r  => "mgn",
         S  => "WSN",  s  => "wsn",
         T  => "ACN",  t  => "acn",
         U  => "TGA",  u  => "tga",
         V  => "GTN",  v  => "gtn",
         W  => "TGG",  w  => "tgg",
         X  => "NNN",  x  => "nnn",
         Y  => "TAY",  y  => "tay",
        '*' => "TRR"
);

%reverse_genetic_code_RNA = (
         A  => "GCN",  a  => "gcn",
         C  => "UGY",  c  => "ugy",
         D  => "GAY",  d  => "gay",
         E  => "GAR",  e  => "gar",
         F  => "UUY",  f  => "uuy",
         G  => "GGN",  g  => "ggn",
         H  => "CAY",  h  => "cay",
         I  => "AUH",  i  => "auh",
         K  => "AAR",  k  => "aar",
         L  => "YUN",  l  => "yun",
         M  => "AUG",  m  => "aug",
         N  => "AAY",  n  => "aay",
         P  => "CCN",  p  => "ccn",
         Q  => "CAR",  q  => "car",
         R  => "MGN",  r  => "mgn",
         S  => "WSN",  s  => "wsn",
         T  => "ACN",  t  => "acn",
         U  => "UGA",  u  => "uga",
         V  => "GUN",  v  => "gun",
         W  => "UGG",  w  => "ugg",
         X  => "NNN",  x  => "nnn",
         Y  => "UAY",  y  => "uay",
        '*' => "URR"
);


%DNA_letter_can_be = (
    A => ["A"],                 a => ["a"],
    B => ["C", "G", "T"],       b => ["c", "g", "t"],
    C => ["C"],                 c => ["c"],
    D => ["A", "G", "T"],       d => ["a", "g", "t"],
    G => ["G"],                 g => ["g"],
    H => ["A", "C", "T"],       h => ["a", "c", "t"],
    K => ["G", "T"],            k => ["g", "t"],
    M => ["A", "C"],            m => ["a", "c"],
    N => ["A", "C", "G", "T"],  n => ["a", "c", "g", "t"],
    R => ["A", "G"],            r => ["a", "g"],
    S => ["C", "G"],            s => ["c", "g"],
    T => ["T"],                 t => ["t"],
    U => ["T"],                 u => ["t"],
    V => ["A", "C", "G"],       v => ["a", "c", "g"],
    W => ["A", "T"],            w => ["a", "t"],
    Y => ["C", "T"],            y => ["c", "t"]
);


%RNA_letter_can_be = (
    A => ["A"],                 a => ["a"],
    B => ["C", "G", "U"],       b => ["c", "g", "u"],
    C => ["C"],                 c => ["c"],
    D => ["A", "G", "U"],       d => ["a", "g", "u"],
    G => ["G"],                 g => ["g"],
    H => ["A", "C", "U"],       h => ["a", "c", "u"],
    K => ["G", "U"],            k => ["g", "u"],
    M => ["A", "C"],            m => ["a", "c"],
    N => ["A", "C", "G", "U"],  n => ["a", "c", "g", "u"],
    R => ["A", "G"],            r => ["a", "g"],
    S => ["C", "G"],            s => ["c", "g"],
    T => ["U"],                 t => ["u"],
    U => ["U"],                 u => ["u"],
    V => ["A", "C", "G"],       v => ["a", "c", "g"],
    W => ["A", "U"],            w => ["a", "u"],
    Y => ["C", "U"],            y => ["c", "u"]
);


%one_letter_to_three_letter_aa = (
         A  => "Ala", a  => "Ala",
         B  => "Asx", b  => "Asx",
         C  => "Cys", c  => "Cys",
         D  => "Asp", d  => "Asp",
         E  => "Glu", e  => "Glu",
         F  => "Phe", f  => "Phe",
         G  => "Gly", g  => "Gly",
         H  => "His", h  => "His",
         I  => "Ile", i  => "Ile",
         K  => "Lys", k  => "Lys",
         L  => "Leu", l  => "Leu",
         M  => "Met", m  => "Met",
         N  => "Asn", n  => "Asn",
         P  => "Pro", p  => "Pro",
         Q  => "Gln", q  => "Gln",
         R  => "Arg", r  => "Arg",
         S  => "Ser", s  => "Ser",
         T  => "Thr", t  => "Thr",
         U  => "Sec", u  => "Sec",
         V  => "Val", v  => "Val",
         W  => "Trp", w  => "Trp",
         X  => "Xxx", x  => "Xxx",
         Y  => "Tyr", y  => "Tyr",
         Z  => "Glx", z  => "Glx",
        '*' => "***"
        );


%three_letter_to_one_letter_aa = (
     ALA  => "A",   Ala  => "A",   ala  => "a",
     ARG  => "R",   Arg  => "R",   arg  => "r",
     ASN  => "N",   Asn  => "N",   asn  => "n",
     ASP  => "D",   Asp  => "D",   asp  => "d",
     ASX  => "B",   Asx  => "B",   asx  => "b",
     CYS  => "C",   Cys  => "C",   cys  => "c",
     GLN  => "Q",   Gln  => "Q",   gln  => "q",
     GLU  => "E",   Glu  => "E",   glu  => "e",
     GLX  => "Z",   Glx  => "Z",   glx  => "z",
     GLY  => "G",   Gly  => "G",   gly  => "g",
     HIS  => "H",   His  => "H",   his  => "h",
     ILE  => "I",   Ile  => "I",   ile  => "i",
     LEU  => "L",   Leu  => "L",   leu  => "l",
     LYS  => "K",   Lys  => "K",   lys  => "k",
     MET  => "M",   Met  => "M",   met  => "m",
     PHE  => "F",   Phe  => "F",   phe  => "f",
     PRO  => "P",   Pro  => "P",   pro  => "p",
     SEC  => "U",   Sec  => "U",   sec  => "u",
     SER  => "S",   Ser  => "S",   ser  => "s",
     THR  => "T",   Thr  => "T",   thr  => "t",
     TRP  => "W",   Trp  => "W",   trp  => "w",
     TYR  => "Y",   Tyr  => "Y",   tyr  => "y",
     VAL  => "V",   Val  => "V",   val  => "v",
     XAA  => "X",   Xaa  => "X",   xaa  => "x",
     XXX  => "X",   Xxx  => "X",   xxx  => "x",
    '***' => "*"
);


#-----------------------------------------------------------------------------
#  Translate nucleotides to one letter protein.  Respects case of the
#  nucleotide sequence.
#
#      $aa = translate_seq( $nt, $met_start )
#      $aa = translate_seq( $nt )
#
#-----------------------------------------------------------------------------

sub translate_seq
{
    my $seq = shift;
    $seq =~ tr/-//d;        #  remove gaps

    my @codons = $seq =~ m/(...?)/g;  #  Will try to translate last 2 nt

    #  A second argument that is true forces first amino acid to be Met

    my @met;
    if ( ( shift @_ ) && ( my $codon1 = shift @codons ) )
    {
        push @met, ( $codon1 =~ /[a-z]/ ? 'm' : 'M' );
    }

    join( '', @met, map { translate_codon( $_ ) } @codons )
}


#-----------------------------------------------------------------------------
#  Translate a single triplet with "universal" genetic code.
#
#      $aa = translate_codon( $triplet )
#      $aa = translate_DNA_codon( $triplet )
#      $aa = translate_uc_DNA_codon( $triplet )
#
#-----------------------------------------------------------------------------

sub translate_DNA_codon { translate_codon( @_ ) }

sub translate_uc_DNA_codon { translate_codon( uc $_[0] ) }

sub translate_codon
{
    my $codon = shift;
    $codon =~ tr/Uu/Tt/;     #  Make it DNA

    #  Try a simple lookup:

    my $aa;
    if ( $aa = $genetic_code{ $codon } ) { return $aa }

    #  Attempt to recover from mixed-case codons:

    $codon = ( $codon =~ /[a-z]/ ) ? lc $codon : uc $codon;
    if ( $aa = $genetic_code{ $codon } ) { return $aa }

    #  The code defined above catches simple N, R and Y ambiguities in the
    #  third position.  Other codons (e.g., GG[KMSWBDHV], or even GG) might
    #  be unambiguously translated by converting the last position to N and
    #  seeing if this is in the code table:

    my $N = ( $codon =~ /[a-z]/ ) ? 'n' : 'N';
    if ( $aa = $genetic_code{ substr($codon,0,2) . $N } ) { return $aa }

    #  Test that codon is valid for an unambiguous aa:

    my $X = ( $codon =~ /[a-z]/ ) ? 'x' : 'X';
    if ( $codon !~ m/^[ACGTMY][ACGT][ACGTKMRSWYBDHVN]$/i
      && $codon !~ m/^YT[AGR]$/i     #  Leu YTR
      && $codon !~ m/^MG[AGR]$/i     #  Arg MGR
       )
    {
        return $X;
    }

    #  Expand all ambiguous nucleotides to see if they all yield same aa.

    my @n1 = @{ $DNA_letter_can_be{ substr( $codon, 0, 1 ) } };
    my $n2 =                        substr( $codon, 1, 1 );
    my @n3 = @{ $DNA_letter_can_be{ substr( $codon, 2, 1 ) } };
    my @triples = map { my $n12 = $_ . $n2; map { $n12 . $_ } @n3 } @n1;

    my $triple = shift @triples;
    $aa = $genetic_code{ $triple };
    $aa or return $X;

    foreach $triple ( @triples ) { return $X if $aa ne $genetic_code{$triple} }

    return $aa;
}


#-----------------------------------------------------------------------------
#  Translate with a user-supplied genetic code to translate a sequence.
#  Diagnose the use of upper versus lower, and T versus U in the supplied
#  code, and transform the supplied nucleotide sequence to match.
#
#     $aa = translate_seq_with_user_code( $nt, \%gen_code )
#     $aa = translate_seq_with_user_code( $nt, \%gen_code, $start_with_met )
#
#  Modified 2007-11-22 to be less intrusive in these diagnoses by sensing
#  the presence of both versions in the user code.
#-----------------------------------------------------------------------------

sub translate_seq_with_user_code
{
    my $seq = shift;
    $seq =~ tr/-//d;     #  remove gaps  ***  Why?

    my $gc = shift;      #  Reference to hash of code
    if (! $gc || ref($gc) ne "HASH")
    {
        print STDERR "translate_seq_with_user_code needs genetic code hash as second argument.";
        return undef;
    }

    #  Test code support for upper vs lower case:

    my ( $TTT, $UUU );
    if    ( $gc->{AAA} && ! $gc->{aaa} )   #  Uppercase only code table
    {
        $seq = uc $seq;     #  Uppercase sequence
        ( $TTT, $UUU ) = ( 'TTT', 'UUU' );
    }
    elsif ( $gc->{aaa} && ! $gc->{AAA} )   #  Lowercase only code table
    {
        $seq = lc $seq;     #  Lowercase sequence
        ( $TTT, $UUU ) = ( 'ttt', 'uuu' );
    }
    elsif ( $gc->{aaa} )
    {
        ( $TTT, $UUU ) = ( 'ttt', 'uuu' );
    }
    else
    {
        print STDERR "User-supplied genetic code does not have aaa or AAA\n";
        return undef;
    }

    #  Test code support for U vs T:

    my $ambigs;
    if    ( $gc->{$UUU} && ! $gc->{$TTT} )  # RNA only code table
    {
        $seq = tr/Tt/Uu/;
        $ambigs = \%RNA_letter_can_be;
    }
    elsif ( $gc->{$TTT} && ! $gc->{$UUU} )  # DNA only code table
    {
        $seq = tr/Uu/Tt/;
        $ambigs = \%DNA_letter_can_be;
    }
    else
    {
        my $t = $seq =~ tr/Tt//;
        my $u = $seq =~ tr/Uu//;
        $ambigs = ( $t > $u ) ? \%DNA_letter_can_be : \%RNA_letter_can_be;
    }

    #  We can now do the codon-by-codon translation:

    my @codons = $seq =~ m/(...?)/g;  #  will try to translate last 2 nt

    #  A third argument that is true forces first amino acid to be Met

    my @met;
    if ( ( shift @_ ) && ( my $codon1 = shift @codons ) )
    {
        push @met, ( $codon1 =~ /[a-z]/ ? 'm' : 'M' );
    }

    join( '', @met, map { translate_codon_with_user_code( $_, $gc, $ambigs ) } @codons )
}


#-----------------------------------------------------------------------------
#  Translate with user-supplied genetic code hash.  No error check on the code.
#  Should only be called through translate_seq_with_user_code.
#
#     $aa = translate_codon_with_user_code( $triplet, \%code, \%ambig_table )
#
#   $triplet      speaks for itself
#  \%code         ref to the hash with the codon translations
#  \%ambig_table  ref to hash with lists of nucleotides for each ambig code
#-----------------------------------------------------------------------------

sub translate_codon_with_user_code
{
    my ( $codon, $gc, $ambigs ) = @_;

    #  Try a simple lookup:

    my $aa;
    if ( $aa = $gc->{ $codon } ) { return $aa }

    #  Attempt to recover from mixed-case codons:

    if ( $aa = $gc->{ lc $codon } ) { return( $gc->{ $codon } = $aa ) }
    if ( $aa = $gc->{ uc $codon } ) { return( $gc->{ $codon } = $aa ) }

    $codon = ( $codon =~ /[a-z]/ ) ? lc $codon : uc $codon;

    #  Test that codon is valid for an unambiguous aa:

    my $X = ( $codon =~ /[a-z]/ ) ? 'x' : 'X';

    if ( $codon =~ m/^[ACGTU][ACGTU]$/i )  # Add N?
    {
        $codon .= ( $codon =~ /[a-z]/ ) ? 'n' : 'N';
    }
    #  This makes assumptions about the user code, but tranlating ambiguous
    #  codons is really a bit off the wall to start with:
    elsif ( $codon !~ m/^[ACGTUMY][ACGTU][ACGTUKMRSWYBDHVN]$/i ) # Valid?
    {
        return( $gc->{ $codon } = $X );
    }

    #  Expand all ambiguous nucleotides to see if they all yield same aa.

    my @n1 = @{ $ambigs->{ substr( $codon, 0, 1 ) } };
    my $n2 =               substr( $codon, 1, 1 );
    my @n3 = @{ $ambigs->{ substr( $codon, 2, 1 ) } };
    my @triples = map { my $n12 = $_ . $n2; map { $n12 . $_ } @n3 } @n1;

    my $triple = shift @triples;
    $aa = $gc->{ $triple } || $gc->{ lc $triple } || $gc->{ uc $triple };
    $aa or return( $gc->{ $codon } = $X );

    foreach $triple ( @triples )
    {
        return( $gc->{ $codon } = $X ) if $aa ne ( $gc->{$triple} || $gc->{lc $triple} || $gc->{uc $triple} );
    }

    return( $gc->{ $codon } = $aa );
}


#-----------------------------------------------------------------------------
#  Read a list of intervals from a file.
#  Allow id_start_end, or id \s start \s end formats
#
#     @intervals = read_intervals( \*FILEHANDLE )
#-----------------------------------------------------------------------------
sub read_intervals {
    my $fh = shift;
    my @intervals = ();

    while (<$fh>) {
        chomp;
           /^(\S+)_(\d+)_(\d+)(\s.*)?$/        #  id_start_end       WIT2
        || /^(\S+)_(\d+)-(\d+)(\s.*)?$/        #  id_start-end       ???
        || /^(\S+)=(\d+)=(\d+)(\s.*)?$/        #  id=start=end       Badger
        || /^(\S+)\s+(\d+)\s+(\d+)(\s.*)?$/    #  id \s start \s end
        || next;

        #  Matched a pattern.  Store reference to (id, left, right):
        push @intervals, ($2 < $3) ? [ $1, $2+0, $3+0 ]
                                   : [ $1, $3+0, $2+0 ];
    }
    return @intervals;
}


#-----------------------------------------------------------------------------
#  Convert a list of intervals to read [ id, left_end, right_end ].
#
#     @intervals = standardize_intervals( @interval_refs )
#-----------------------------------------------------------------------------
sub standardize_intervals {
    map { ( $_->[1] < $_->[2] ) ? $_ : [ $_->[0], $_->[2], $_->[1] ] } @_;
}


#-----------------------------------------------------------------------------
#  Take the union of a list of intervals
#
#     @joined = join_intervals( @interval_refs )
#-----------------------------------------------------------------------------
sub join_intervals {
    my @ordered = sort { $a->[0] cmp $b->[0]   # first by id
                      || $a->[1] <=> $b->[1]   # next by left end
                      || $b->[2] <=> $a->[2]   # finally longest first
                       } @_;

    my @joined = ();
    my $n_int = @ordered;

    my ($cur_id)    = "";
    my ($cur_left)  = -1;
    my ($cur_right) = -1;
    my ($new_id, $new_left, $new_right);

    for (my $i = 0; $i < $n_int; $i++) {
        ($new_id, $new_left, $new_right) = @{$ordered[$i]};  # get the new data

        if ( ( $new_id ne $cur_id)          # if new contig
          || ( $new_left > $cur_right + 1)  #    or not touching previous
           ) {                              # push the previous interval
            if ($cur_id) { push (@joined, [ $cur_id, $cur_left, $cur_right ]) }
            $cur_id = $new_id;              # update the current interval
            $cur_left = $new_left;
            $cur_right = $new_right;
        }

        elsif ($new_right > $cur_right) {   # extend the right end if necessary
            $cur_right = $new_right;
        }
    }

    if ($cur_id) { push (@joined, [$cur_id, $cur_left, $cur_right]) }
    return @joined;
}


#-----------------------------------------------------------------------------
#  Split location strings to oriented intervals.
#
#     @intervals = locations_2_intervals( @locations )
#     $interval  = locations_2_intervals( $location  )
#-----------------------------------------------------------------------------
sub locations_2_intervals {
    my @intervals = map { /^(\S+)_(\d+)_(\d+)(\s.*)?$/
                       || /^(\S+)_(\d+)-(\d+)(\s.*)?$/
                       || /^(\S+)=(\d+)=(\d+)(\s.*)?$/
                       || /^(\S+)\s+(\d+)\s+(\d+)(\s.*)?$/
                        ? [ $1, $2+0, $3+0 ]
                        : ()
                        } @_;

    return wantarray ? @intervals : $intervals[0];
}


#-----------------------------------------------------------------------------
#  Read a list of oriented intervals from a file.
#  Allow id_start_end, or id \s start \s end formats
#
#     @intervals = read_oriented_intervals( \*FILEHANDLE )
#-----------------------------------------------------------------------------
sub read_oriented_intervals {
    my $fh = shift;
    my @intervals = ();

    while (<$fh>) {
        chomp;
           /^(\S+)_(\d+)_(\d+)(\s.*)?$/        #  id_start_end       WIT2
        || /^(\S+)_(\d+)-(\d+)(\s.*)?$/        #  id_start-end       ???
        || /^(\S+)=(\d+)=(\d+)(\s.*)?$/        #  id=start=end       Badger
        || /^(\S+)\s+(\d+)\s+(\d+)(\s.*)?$/    #  id \s start \s end
        || next;

        #  Matched a pattern.  Store reference to (id, start, end):
        push @intervals, [ $1, $2+0, $3+0 ];
    }
    return @intervals;
}


#-----------------------------------------------------------------------------
#  Reverse the orientation of a list of intervals
#
#     @reversed = reverse_intervals( @interval_refs )
#-----------------------------------------------------------------------------
sub reverse_intervals {
    map { [ $_->[0], $_->[2], $_->[1] ] } @_;
}


#-----------------------------------------------------------------------------
#  Convert GenBank locations to SEED locations
#
#     @seed_locs = gb_location_2_seed( $contig, @gb_locs )
#-----------------------------------------------------------------------------
sub gb_location_2_seed
{
    my $contig = shift @_;
    $contig or die "First arg of gb_location_2_seed must be contig_id\n";

    map { join( ',', gb_loc_2_seed_2( $contig, $_ ) ) || undef } @_
}

sub gb_loc_2_seed_2
{
    my ( $contig, $loc ) = @_;

    if ( $loc =~ /^(\d+)\.\.(\d+)$/ )
    {
        join( '_', $contig, $1, $2 )
    }

    elsif ( $loc =~ /^join\((.*)\)$/ )
    {
        $loc = $1;
        my $lvl = 0;
        for ( my $i = length( $loc )-1; $i >= 0; $i-- )
        {
            for ( substr( $loc, $i, 1 ) )
            {
                /,/ && ! $lvl and substr( $loc, $i, 1 ) = "\t";
                /\(/          and $lvl--;
                /\)/          and $lvl++;
            }
        }
        $lvl == 0 or print STDERR "Paren matching error: $loc\n" and die;
        map { gb_loc_2_seed_2( $contig, $_ ) } split /\t/, $loc
    }

    elsif ( $loc =~ /^complement\((.*)\)$/ )
    {
        map { s/_(\d+)_(\d+)$/_$2_$1/; $_ }
        reverse
        gb_loc_2_seed_2( $contig, $1 )
    }

    else
    {
        ()
    }
}


#-----------------------------------------------------------------------------
#  Read qual.
#
#  Save the contents in a list of refs to arrays: [ $id, $descript, \@qual ]
#
#     @seq_entries = read_qual( )               #  STDIN
#    \@seq_entries = read_qual( )               #  STDIN
#     @seq_entries = read_qual( \*FILEHANDLE )
#    \@seq_entries = read_qual( \*FILEHANDLE )
#     @seq_entries = read_qual(  $filename )
#    \@seq_entries = read_qual(  $filename )
#-----------------------------------------------------------------------------
sub read_qual {
    my ( $fh, $name, $close, $unused ) = input_filehandle( $_[0] );
    $unused && die "Bad reference type (" . ref( $unused ) . ") passed to read_qual\n";

    my @quals = ();
    my ($id, $desc, $qual) = ("", "", []);

    while ( <$fh> ) {
        chomp;
        if (/^>\s*(\S+)(\s+(.*))?$/) {        #  new id
            if ($id && @$qual) { push @quals, [ $id, $desc, $qual ] }
            ($id, $desc, $qual) = ($1, $3 ? $3 : "", []);
        }
        else {
            push @$qual, split;
        }
    }
    close( $fh ) if $close;

    if ($id && @$qual) { push @quals, [ $id, $desc, $qual ] }
    return wantarray ? @quals : \@quals;
}


#-------------------------------------------------------------------------------
#  Fraction difference for an alignment of two nucleotide sequences in terms of
#  number of differing residues, number of gaps, and number of gap opennings.
#
#     $fraction_diff = fraction_nt_diff( $seq1, $seq2, \%options )
#
#  or
#
#     $fraction_diff = fraction_nt_diff( $seq1, $seq2 )
#     $fraction_diff = fraction_nt_diff( $seq1, $seq2, $gap_wgt )
#     $fraction_diff = fraction_nt_diff( $seq1, $seq2, $open_wgt, $extend_wgt )
#
#  Options:
#
#      gap      => $gap_wgt          # Gap open and extend weight (D = 0.5)
#      open     => $open_wgt         # Gap openning weight (D = gap_wgt)
#      extend   => $extend_wgt       # Gap extension weight (D = open_wgt)
#      t_gap    => $term_gap_wgt     # Terminal open and extend weight
#      t_open   => $term_open_wgt    # Terminal gap open weight (D = open_wgt)
#      t_extend => $term_extend_wgt  # Terminal gap extend weight (D = extend_wgt)
#
#  Default gap open and gap extend weights are 1/2.  Beware that
#
#     $fraction_diff = fraction_nt_diff( $seq1, $seq2, 1 )
#
#  and
#
#     $fraction_diff = fraction_nt_diff( $seq1, $seq2, 1, 0 )
#
#  are different.  The first has equal openning and extension weights, whereas
#  the second has an openning weight of 1, and and extension weight of 0 (it
#  only penalizes the number of runs of gaps).
#-------------------------------------------------------------------------------
sub fraction_nt_diff
{
    my ( $npos, $nid, $ndif, $ngap, $nopen, $tgap, $topen ) = interpret_nt_align( @_[0,1] );

    my $diff_scr;
    if ( ref( $_[2] ) eq 'HASH' )
    {
        my $opts = $_[2];
        my $gap_open    = defined $opts->{ open }     ? $opts->{ open }
                        : defined $opts->{ gap }      ? $opts->{ gap }
                        : 0.5;
        my $gap_extend  = defined $opts->{ extend }   ? $opts->{ extend }
                        : $gap_open;
        my $term_open   = defined $opts->{ t_open }   ? $opts->{ t_open }
                        : defined $opts->{ t_gap }    ? $opts->{ t_gap }
                        : $gap_open;
        my $term_extend = defined $opts->{ t_extend } ? $opts->{ t_extend }
                        : defined $opts->{ t_gap }    ? $opts->{ t_gap }
                        : $gap_extend;

        $nopen -= $topen;
        $ngap  -= $tgap;
        $diff_scr = $ndif + $gap_open  * $nopen + $gap_extend  * ($ngap-$nopen)
                          + $term_open * $topen + $term_extend * ($tgap-$topen);
    }
    else
    {
        my $gap_open   = defined( $_[2] ) ? $_[2] : 0.5;
        my $gap_extend = defined( $_[3] ) ? $_[3] : $gap_open;
        $diff_scr = $ndif + $gap_open * $nopen + $gap_extend * ($ngap-$nopen);
    }
    my $ttl_scr = $nid + $diff_scr;

    $ttl_scr ? $diff_scr / $ttl_scr : undef
}


#-------------------------------------------------------------------------------
#  Interpret an alignment of two nucleotide sequences in terms of: useful
#  aligned positions (unambiguous, and not a common gap), number of identical
#  residues, number of differing residues, number of gaps, and number of gap
#  opennings.
#
#     ( $npos, $nid, $ndif, $ngap, $nopen, $tgap, $topen ) = interpret_nt_align( $seq1, $seq2 )
#
#  $npos  = total aligned positons (= $nid + $ndif + $ngap)
#  $nid   = number of positions with identical nucleotides (ignoring case)
#  $ndif  = number of positions with differing nucleotides
#  $ngap  = number of positions with gap in one sequence but not the other
#  $nopen = number of runs of gaps
#  $tgap  = number of gaps in runs adjacent to a terminus
#  $topen = number of alignment ends with gaps
#
#  Some of the methods might seem overly complex, but are necessary for cases
#  in which the gaps switch strands in the alignment:
#
#     seq1  ---ACGTGAC--TTGCAGAG
#     seq2  TTT---TGACGG--GCAGGG
#     mask  00000011110000111111
#
#     npos  = 20
#     nid   =  9
#     ndif  =  1
#     ngap  = 10
#     nopen =  4
#     tgap  =  3
#     topen =  1
#
#  Although there are 4 gap opennings, there are only 2 runs in the mask,
#  and the terminal run is length 6, not 3.  (Why handle these?  Because
#  pairs of sequences from a multiple sequence alignment can look like this.)
#-------------------------------------------------------------------------------
sub interpret_nt_align
{
    #  Remove alignment columns that are not informative:
    my ( $s1, $s2 ) = useful_nt_align( @_[0,1] );
    my $nmat = length( $s1 );          # Useful alignment length

    my $m1 = $s1;
    $m1 =~ tr/ACGT/\377/;              # Nucleotides to all 1 bits
    $m1 =~ tr/\377/\000/c;             # Others (gaps) to null byte
    my $m2 = $s2;
    $m2 =~ tr/ACGT/\377/;              # Nucleotides to all 1 bits
    $m2 =~ tr/\377/\000/c;             # Others (gaps) to null byte
    $m1 &= $m2;                        # Gap in either sequence becomes null
    $s1 &= $m1;                        # Apply mask to sequence 1
    $s2 &= $m1;                        # Apply mask to sequence 2
    my $nopen = @{[ $s1 =~ /\000+/g ]}   # Gap opens in sequence 1
              + @{[ $s2 =~ /\000+/g ]};  # Gap opens in sequence 2
    my ( $tgap, $topen ) = ( 0, 0 );
    if ( $s1 =~ /^(\000+)/ || $s2 =~ /^(\000+)/ ) { $tgap += length( $1 ); $topen++ }
    if ( $s1 =~ /(\000+)$/ || $s2 =~ /(\000+)$/ ) { $tgap += length( $1 ); $topen++ }
    $s1 =~ tr/\000//d;                 # Remove nulls (former gaps)
    $s2 =~ tr/\000//d;                 # Remove nulls (former gaps)
    my $ngap = $nmat - length( $s1 );  # Total gaps

    my $xor = $s1 ^ $s2;               # xor of identical residues is null byte
    my $nid = ( $xor =~ tr/\000//d );  # Count the nulls (identical residues)
    my $ndif = $nmat - $nid - $ngap;

    ( $nmat, $nid, $ndif, $ngap, $nopen, $tgap, $topen )
}


sub useful_nt_align
{
    my ( $s1, $s2 ) = map { uc $_ } @_;
    $s1 =~ tr/U/T/;         # Convert U to T
    my $m1 = $s1;
    $m1 =~ tr/ACGT-/\377/;  # Allowed symbols to hex FF byte
    $m1 =~ tr/\377/\000/c;  # All else to null byte
    $s2 =~ tr/U/T/;         # Convert U to T
    my $m2 = $s2;
    $m2 =~ tr/ACGT-/\377/;  # Allowed symbols to hex FF byte
    $m2 =~ tr/\377/\000/c;  # All else to null byte
    $m1 &= $m2;             # Invalid in either sequence becomes null
    $s1 &= $m1;             # Apply mask to sequence 1
    $s1 =~ tr/\000//d;      # Delete nulls in sequence 1
    $s2 &= $m1;             # Apply mask to sequence 2
    $s2 =~ tr/\000//d;      # Delete nulls in sequence 2
    ( $s1, $s2 )
}


#-------------------------------------------------------------------------------
#  Interpret an alignment of two protein sequences in terms of: useful
#  aligned positions (unambiguous, and not a common gap), number of identical
#  residues, number of differing residues, number of gaps, and number of gap
#  opennings.
#
#     ( $npos, $nid, $ndif, $ngap, $nopen, $tgap, $topen ) = interpret_aa_align( $seq1, $seq2 )
#
#  $npos  = total aligned positons (= $nid + $ndif + $ngap)
#  $nid   = number of positions with identical amino acids (ignoring case)
#  $ndif  = number of positions with differing amino acids
#  $ngap  = number of positions with gap in one sequence but not the other
#  $nopen = number of runs of gaps
#  $tgap  = number of gaps in runs adjacent to a terminus
#  $topen = number of alignment ends with gaps
#
#-------------------------------------------------------------------------------
sub interpret_aa_align
{
    #  Remove alignment columns that are not informative:
    my ( $s1, $s2 ) = useful_aa_align( @_[0,1] );
    my $nmat = length( $s1 );            # Useful alignment length

    my $m1 = $s1;
    $m1 =~ tr/ACDEFGHIKLMNPQRSTUVWY/\377/;  # Amino acids to all 1 bits
    $m1 =~ tr/\377/\000/c;               # Others (gaps) to null byte
    my $m2 = $s2;
    $m2 =~ tr/ACDEFGHIKLMNPQRSTUVWY/\377/;  # Amino acids to all 1 bits
    $m2 =~ tr/\377/\000/c;               # Others (gaps) to null byte
    $m1 &= $m2;                          # Gap in either sequence becomes null
    $s1 &= $m1;                          # Apply mask to sequence 1
    $s2 &= $m1;                          # Apply mask to sequence 2
    my $nopen = @{[ $s1 =~ /\000+/g ]}   # Gap opens in sequence 1
              + @{[ $s2 =~ /\000+/g ]};  # Gap opens in sequence 2
    my ( $tgap, $topen ) = ( 0, 0 );
    if ( $s1 =~ /^(\000+)/ || $s2 =~ /^(\000+)/ ) { $tgap += length( $1 ); $topen++ }
    if ( $s1 =~ /(\000+)$/ || $s2 =~ /(\000+)$/ ) { $tgap += length( $1 ); $topen++ }
    $s1 =~ tr/\000//d;                 # Remove nulls (former gaps)
    $s2 =~ tr/\000//d;                 # Remove nulls (former gaps)
    my $ngap = $nmat - length( $s1 );  # Total gaps

    my $xor = $s1 ^ $s2;               # xor of identical residues is null byte
    my $nid = ( $xor =~ tr/\000//d );  # Count the nulls (identical residues)
    my $ndif = $nmat - $nid - $ngap;

    ( $nmat, $nid, $ndif, $ngap, $nopen, $tgap, $topen )
}


sub useful_aa_align
{
    my ( $s1, $s2 ) = map { uc $_ } @_;
    my $m1 = $s1;
    $m1 =~ tr/ACDEFGHIKLMNPQRSTUVWY-/\377/;  # Allowed symbols to hex FF byte
    $m1 =~ tr/\377/\000/c;  # All else to null byte
    my $m2 = $s2;
    $m2 =~ tr/ACDEFGHIKLMNPQRSTUVWY-/\377/;  # Allowed symbols to hex FF byte
    $m2 =~ tr/\377/\000/c;  # All else to null byte
    $m1 &= $m2;             # Invalid in either sequence becomes null
    $s1 &= $m1;             # Apply mask to sequence 1
    $s1 =~ tr/\000//d;      # Delete nulls in sequence 1
    $s2 &= $m1;             # Apply mask to sequence 2
    $s2 =~ tr/\000//d;      # Delete nulls in sequence 2
    ( $s1, $s2 )
}


#-------------------------------------------------------------------------------
#  Return the fraction identity for oligomers over a range of lengths
#
#     @sims = oligomer_similarity( $seq1, $seq2, \%opts )
#
#-------------------------------------------------------------------------------
sub oligomer_similarity
{
    my ( $seq1, $seq2, $opts ) = @_;
    $seq1 && $seq2 or return ();
    $seq1 = $seq1->[2] if ref( $seq1 ) eq 'ARRAY';
    $seq2 = $seq2->[2] if ref( $seq2 ) eq 'ARRAY';
    $opts && ref( $opts ) eq 'HASH' or ( $opts = {} );

    my $min = $opts->{ min } || $opts->{ max } || 2;
    my $max = $opts->{ max } || $min;
    return () if $max < $min;

    #  Remove shared gaps

    my $mask1 = gap_mask( $seq1 );
    my $mask2 = gap_mask( $seq2 );
    my $mask  = $mask1 | $mask2;
    $seq1     = $seq1 & $mask;          # Apply mask to sequence
    $seq1     =~ tr/\000//d;     # Delete null characters
    $seq2     = $seq2 & $mask;          # Apply mask to sequence
    $seq2     =~ tr/\000//d;     # Delete null characters
    #  Remove terminal gaps

    $mask  = $mask1 & $mask2;
    my $n1 = $mask =~ /^(\000+)/ ? length( $1 ) : 0;
    my $n2 = $mask =~ /(\000+)$/ ? length( $1 ) : 0;
    if ( $n1 || $n2 )
    {
        my $len = length( $seq1 ) - ( $n1 + $n2 );
        $seq1 = substr( $seq1, $n1, $len );
        $seq2 = substr( $seq2, $n1, $len );
    }

    #  Find the runs of identity

    my $xor = $seq1 ^ $seq2;           # xor of identical residues is null byte

    # Runs with one or more identities

    my %cnt;
    foreach ( $xor =~ m/(\000+)/g ) { $cnt{ length($_) }++ }
    map { my $n = $_;
          my $ttl = 0;
          for ( grep { $_ >= $n } keys %cnt ) { $ttl += $cnt{$_} * ( $_ - ($n-1) ) }
          my $nmax = length( $xor ) - ($n-1); 
          $nmax > 0 ? $ttl / $nmax : undef;
        } ( $min .. $max );
}


1;
