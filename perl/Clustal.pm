#__perl__

use strict;
package Clustal;
use Data::Dumper;

=pod

=head1 Clustal

A lightweight clustal module written by Rob Edwards. Really only reads a clustal format alignment and gets those sequences back

=cut

=head1 new

Instantiate a new method. You can pass in a file
my $clustal = Clustal->new(-file => "out.aln");

=cut

sub new {
	my ($class, %args)=@_;
	my $self         = {};
	if ($args{'-file'} && -e $args{'-file'}) {
		$self->{file} = $args{'-file'};
	}

	return bless $self, $class;
}


=head1 file

Set the alignment file

$clustal->file($file);

=cut

sub file {
	my ($self, $file) = @_;
	if ($file) {
		if (-e $file) {
			$self->{file} = $file;
		} else {
			die "$file not found";
		}
	}
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	return $self->{file};
}

=head1 version

Get the clustal version used in the alignment

$ver = $clustal->version()

=cut 

sub version {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	if (!defined $self->{version}) {$self->parse_file()}
	return $self->{version};
}


=head1 parse_file

Parse a file. You don't need to call this if you call the alignment() method.

=cut

sub parse_file {
	my ($self)=@_;
	unless ($self->{'file'}) {die "No clustal file is specified. Please set a file before you try and parse it"}
	return if ($self->{'parsed_file'} && $self->{'parsed_file'} eq $self->{'file'}); # nothing to do!
	$self->{'parsed_file'} = $self->{'file'};
	open(IN, $self->{'parsed_file'}) || die "Can't open " . $self->{'parsed_file'};
	my $header = <IN>;
	if ($header =~ /CLUSTAL (\S+) multiple sequence alignment/) {$self->{'version'}=$1}
	else {die $self->{'parsed_file'}. " does not appear to be a clustal file"}
	my $currindex=0; my $spaces; 
	$self->{'identities'} = undef;
	$self->{'index'} = undef;
	$self->{'sequences'} = undef;

	while (<IN>) {
		chomp;
		next if (/^\s+$/);
		if (!$spaces && /^(\S+\s+)\S+/) {$spaces = length($1)} # how many spaces are there before the alignment
		if (/^[\s\:\*\.]+/) {
			# this is the identity line
			s/^.{$spaces}//;
			$self->{'identities'} .= $_;
			next;
		}
		if (/^(\S+)\s+(\S+)$/) {
			my ($id, $seq)=($1, $2);
			unless (defined $self->{'index'}->{$id}) {$self->{'index'}->{$id}=$currindex; $currindex++}
			$self->{'sequences'}->{$id}.=$seq;
		}
	}
	close IN;
}


=head1 ids

Get the sequence IDs in the order that they appear in the file.

my @ids = $clustalw->id();

=cut

sub ids {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	if (!defined $self->{'index'}) {$self->parse_file()}
	my @ids = sort {$self->{'index'}->{$a} <=> $self->{'index'}->{$b}} keys %{$self->{'index'}};
	return @ids;
}

=head1 get_a_sequence

Get the sequence for a single ID

my $sequence = $clustal->get_a_sequence($id);

=cut

sub get_a_sequence {
	my ($self, $id) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	if (!defined $self->{'index'}) {$self->parse_file()}
	if (!defined $self->{'index'}->{$id}) {die "$id does not appear in the alignment"}
	return $self->{'sequences'}->{$id};
}

=head1 alignments

Get all the sequences in the file

my @alignment = $clustal->alignments();

This returns an array of arrays. Each element in the array is an array of the sequences.
Thus, $alignment->[0]->[0] aligns with $aligment->[1]->[0]; $alignment->[0]->[1]aligns with $aligment->[1]->[1], etc

=cut

sub alignments {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	my @result;
	foreach my $id ($self->ids()) {
		my @seq = split //, $self->get_a_sequence($id);
		push @result, \@seq;
	}
	return \@result;
}


=head1 trimmed_alignments

The trimmed alignments remove contiguous gaps at the 5' and 3' end of the sequences

my @trimmed_alignment = $clustal->trimmed_alignments()

=cut

sub trimmed_alignments {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	my @result;
	my $fivetrim = 0; my $threetrim = 0;
	# first iterate through and get the lengths to trim
	foreach my $id ($self->ids()) {
		my $seq = $self->get_a_sequence($id);
		if ($seq =~ /^(\-+)/) {
			(length($1) > $fivetrim) ? ($fivetrim = length($1)) : 1;
		}
		if ($seq =~ /(\-+)$/) {
			(length($1) > $threetrim) ? ($threetrim = length($1)) : 1;
		}
	}
	foreach my $id ($self->ids()) {
		my $seq = $self->get_a_sequence($id);
		if ($threetrim) {$seq = substr($seq, 0, length($seq)-$threetrim)}
		if ($fivetrim)  {$seq = substr($seq, $fivetrim)}
		push @result, [split //, $seq];
		$self->{'trimmed_sequence'}->{$id}=$seq;
	}
		
	my $idents = $self->{'identities'};
	if ($threetrim) {$idents = substr($idents, 0, length($idents)-$threetrim)}
	if ($fivetrim && length($idents) > $fivetrim)  {$idents = substr($idents, $fivetrim)}
	$self->{'trimmed_identities'}=$idents;

	return \@result;
}


=head1 identities

Get the identities of the sequences. Returns a reference to an array of all the identities

my $identities = $clustal->identities();

=cut

sub identities {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	if (!defined $self->{'identities'}) {$self->parse_file()}
	return [split //, $self->{'identities'}];
}

=head1 trimmed_identities

Get the percent identity of the 5' and 3' trimmed sequences

my $trimmed_identities = $clustal->trimmed_identities();

=cut

sub trimmed_identities {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	if (!defined $self->{'trimmed_identities'}) {$self->trimmed_alignments()}
	return [split //, $self->{'trimmed_identities'}];
}


=head1 percent_identical

Get the percent of the sequences that are identical. By default this is from the whole sequence, 
but you can also trim off blank sequences at either end, which are then ignored.

my $percent = $clustal->percent_identical()

=cut

sub percent_identical {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	if (!defined $self->{'identities'}) {$self->parse_file()}

	return $self->calc_percent($self->alignments());
}



sub calc_percent {
	my ($self, $alignments)=@_;
	my $n=0; my $same=0;
	for my $j (0 .. $#{$alignments->[0]}) {
		$n++;
		my %bases;
		for my $i (0 .. $#$alignments) {$bases{uc($alignments->[$i]->[$j])}=1}
		(scalar(keys %bases) == 1) && ($same++);
	}
	return sprintf("%.3f", (($same/$n) *100));
}


=head1 trimmed_percent_identical

Get the percents of the trimmed sequences that are identical.

my $percent = $clustal->trimmed_percent_identical()

=cut

sub trimmed_percent_identical {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	if (!defined $self->{'trimmed_identities'}) {$self->trimmed_alignments()}

	return $self->calc_percent($self->trimmed_alignments());
}




=head1 alignment_length

The total length of the alignment, including gaps

my $len = $clustal->alignment_length()

=cut

sub alignment_length {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	my @ids = $self->ids();
	my $seq = $self->get_a_sequence($ids[0]);
	return length($seq);
}

=head1 trimmed_alignment_length

The total length of the trimmed alignment, not including leading and trailing gaps

my $len = $clustal->trimmed_alignment_length()

=cut

sub trimmed_alignment_length {
	my ($self) = @_;
	if ($self->{'parsed_file'} && $self->{'parsed_file'} ne $self->{'file'}) {$self->parse_file()}
	my @ids = $self->ids();
	my $seq = $self->{'trimmed_sequence'}->{$ids[0]};
	return length($seq);
}


1;

