# a collection of simple seed stuff I need to access data from individual directories

package MinSeed;
use strict;
use Rob;
my $rob=new Rob;
use DB_File;
use FileHandle;

=head1 new

Instantiate a min seed. Note that you can give it a directory where all the organisms are expected to be housed:
my $minseed = MinSeed->new({-dir => "/data/Vibrio/rast_annotations/seed_directories/"});

=cut

sub new {
	my ($class, $args)=@_;
	my $self         = {};
	if ($args->{'-dir'} && -d $args->{'-dir'}) {
		$self->{orgdir} = $args->{'-dir'}
	}

	return bless $self, $class;
}

=head1 genomes

An alias for all_genomes

=cut

sub genomes {
	my ($self)=@_;
	return $self->all_genomes();
}

=head1 all_genomes

my @genomes = $minseed->all_genomes();

=cut

sub all_genomes {
	my ($self)=@_;
	return () unless ($self->{orgdir});
	return @{$self->{all_genomes}} if (defined $self->{all_genomes});
	opendir(DIR, $self->{orgdir}) || die "Can't open dir $self->{orgdir}";
	@{$self->{all_genomes}} = grep {-e  $self->{orgdir}."/$_/GENOME"} readdir(DIR);
	return @{$self->{all_genomes}};
}

=head1 is_ion_genome

This is specific for the Vibrio project. I made a file ION in the genome directory that signifies which genomes are from
the Ion Torrent sequencing project.

my @genomes = grep {$minseed->is_ion_genome($_)} $minseed->all_genomes;

if ($minseed->is_ion_genome($genome)) { .. }

=cut 

sub is_ion_genome {
	my ($self, $genome)=@_;

	return undef unless ($genome);
	return undef unless ($self->{orgdir});
	return (-e $self->{orgdir}."/$genome/ION");
}


=head1 genome_of

	get the genome for a given peg
	my $genome = $minseed->genome_of($peg)

=cut

sub genome_of {
	my ($self, $peg) = @_;
	$peg =~ /fig\|(\d+\.\d+)/;
	return $1;
}

=head1 genus_species

Get the genus species of the genome
my $gs = $minseed->genus_species($genome)

=cut

sub genus_species {
	my ($self, $genome)=@_;
	return $self->{$genome}->{genus_species} if (defined $self->{$genome}->{genus_species});
	return unless ($self->{orgdir});
	return unless (-e $self->{orgdir}."/$genome/GENOME");
	open(IN, $self->{orgdir}."/$genome/GENOME") || die "can't open ". $self->{orgdir}."/$genome/GENOME";
	my $g = <IN>;
	chomp $g;
	$self->{$genome}->{genus_species} = $g;
	return $self->{$genome}->{genus_species};
}

=head1 abbreviation

Get the abbreviation for the gneome
my $abbr = $minseed->abbreviation($genome)

=cut

sub abbreviation {
	my ($self, $genome)=@_;
	return $self->{$genome}->{abbreviation} if (defined $self->{$genome}->{abbreviation});
	return unless ($self->{orgdir});
	return $genome unless (-e $self->{orgdir}."/$genome/ABBREVIATION");
	open(IN, $self->{orgdir}."/$genome/ABBREVIATION") || die "can't open ". $self->{orgdir}."/$genome/ABBREVIATION";

	my $g = <IN>;
	chomp $g;
	$self->{$genome}->{abbreviation} = $g;
	return $self->{$genome}->{abbreviation};
}



=head1 pegs_of 

Get the pegs of a genome
my @pegs = $minseed->pegs_of($genome);

=cut

sub pegs_of {
	my ($self, $genome) = @_;
	return keys %{$self->{$genome}->{feature_locations}} if (defined $self->{$genome}->{feature_locations});
	if ($self->{orgdir}) {
		$self->{'last_id'}->{$genome}=0;
		open(TBL, $self->{orgdir}."/$genome/Features/peg/tbl") || die "can't open ".$self->{orgdir}."/$genome/Features/peg/tbl";
		while (<TBL>) {
			chomp;
			my ($peg, $loc) = split /\t/;
			$self->{$genome}->{feature_locations}->{$peg}=$loc;
			$peg =~ /(\d+)$/;
			($1 > $self->{'last_id'}->{$genome}) && ($self->{'last_id'}->{$genome} = $1);
		}
		close TBL;
		# now read the deleted features
		if (-e $self->{orgdir}."/$genome/Features/peg/deleted.features") {
			open(DEL, $self->{orgdir}."/$genome/Features/peg/deleted.features") || die "can't open deleted features: ";
			while (<DEL>) {
				chomp;
				delete $self->{$genome}->{feature_locations}->{$_};
				$self->{$genome}->{is_deleted}->{$_}=1;
			}
			close DEL;
		}

	}
	return keys %{$self->{$genome}->{feature_locations}};
}

=head1 delete_feature

Mark a feature as deleted

=cut

sub delete_feature {
	my ($self, $peg)=@_;
	my $ftype = $self->ftype($peg);
	unless ($ftype) {die "no ftype for $peg"}
	my $genome = $self->genome_of($peg);
	open(DEL, ">>".$self->{orgdir}."/$genome/Features/$ftype/deleted.features") || die "can't open deleted features: $genome/Features/$ftype/deleted.features";
	print DEL $peg, "\n";
	close DEL;
}

=head1 ftype

Figure out the type of a feature (usu rna, peg, etc)

=cut

sub ftype {
	my ($self, $peg) = @_;
	$peg =~ m/fig\|\d+\.\d+\.(\w+)\.\d+/;
	return $1;
}


=head1 is_deleted_fid

Test whether a feature id is deleted

my $test = $fig->is_deleted_fid($peg);

=cut

sub is_deleted_fid {
	my ($self, $fid)=@_;
	my $genome = $self->genome_of($fid);
	if (!defined $self->{$genome}->{feature_locations}) {$self->pegs_of($genome)}
	return $self->{$genome}->{is_deleted}->{$fid}
}

=head1 deleted_fids

Get all deleted fids for a genome

my @deleted_fids = $fig->deleted_fids($genome);

=cut

sub deleted_fids {
	my ($self, $genome)=@_;
	if (!defined $self->{$genome}->{feature_locations}) {$self->pegs_of($genome)}
	return (keys %{$self->{$genome}->{is_deleted}});
}

=head1 feature_location

get the location of a feature

my $loc = $minseed->feature_location($peg);
my @loc = $minseed->feature_location($peg)

In an array context returns each location. In a scalar context returns a comma separated list

=cut

sub feature_location {
	my ($self, $peg) = @_;
	my $genome = $self->genome_of($peg);
	if ($self->{$genome}->{is_deleted}->{$peg}) {print STDERR "WARNING: $peg is deleted\n"}
	unless (defined $self->{$genome}->{feature_locations}) {$self->pegs_of($self->genome_of($peg))}
	unless (defined $self->{$genome}->{feature_locations}->{$peg}) {
		print STDERR "No feature location for $peg ??\n";
		return wantarray ? () : "";
	}
	return wantarray ? split /,/, $self->{$genome}->{feature_locations}->{$peg} : $self->{$genome}->{feature_locations}->{$peg};
}

=head1 set_feature_location

Set a feature location for a given peg (we just append it, then when we read the locations we'll overwrite any existing ones

my $worked = $fig->set_feature_location($peg, $location);

=cut

sub set_feature_location {
	my ($self, $peg, $location) = @_;
	$self->{$self->genome_of($peg)}->{feature_locations}->{$peg}=$location;
	my $type = $self->ftype($peg);
	unless (defined $type) {die "Could not get a type for $peg??"}
	my $genome = $self->genome_of($peg);
	open(TBL, ">>".$self->{orgdir}."/$genome/Features/$type/tbl") || die "can't open ".$self->{orgdir}."/$genome/Features/$type/tbl";
	print TBL "$peg\t$location\n";
	close TBL;
	return 1;
}


=head1 get_dna_seq

get the dna sequence for a protein and additional DNA on either side if required

my $dna = $minseed->get_dna_seq($peg, $extra);



=cut

sub get_dna_seq {
	my ($self, $peg, $extra) = @_;
	my $genome    = $self->genome_of( $peg );
	if ($self->{$genome}->{is_deleted}->{$peg}) {print STDERR "WARNING: $peg is deleted\n"}
	my @locations = $self->feature_location( $peg );
	my $seq = $self->dna_seq($genome, \@locations, $extra);

	return $seq;
}

=head1 dna_seq

usage: $seq = $fig->dna_seq($genome, \@locations, $extra)

Returns the concatenated subsequences described by the reference to the list of locations.  Each location
must be of the form

Contig_Beg_End

where Contig must be the ID of a contig for genome $genome.  If Beg > End the location
describes a stretch of the complementary strand.

If $extra is specified we'll get additional DNA on either side.

=cut

sub dna_seq {
	my($self,$genome,$locations, $extra) = @_;

	my @locations = map { split(/,/,$_) } @$locations;
	my @pieces = ();
	foreach my $loc (@locations)
	{
		if ($loc =~ /^(\S+)_(\d+)_(\d+)$/)
		{
			my ($contig,$beg,$end) = ($1,$2,$3);
			my $ln = $self->contig_ln($genome,$contig);
			if (! $ln) {
				die "dna_seq($genome, $loc): contig length undefined";
			}
			if ($self->between(1,$beg,$ln) && $self->between(1,$end,$ln))
			{
				if ($beg < $end)
				{
					if ($extra) {
						$beg -= $extra;
						($beg < 0) ? $beg = 0 : 1;
						$end += $extra;
						($end > $ln) ? $end += $ln : 1;
					}
					push(@pieces, $self->get_dna($genome,$contig,$beg,$end));
				}
				else
				{
					if ($extra) {
						$end -= $extra;
						($end < 0) ? $end = 0 : 1;
						$beg += $extra;
						($beg > $ln) ? $beg += $ln : 1;
					}
					push(@pieces, $rob->rc($self->get_dna($genome,$contig,$end,$beg)));
				}
			}
		}
	}
	return lc(join("",@pieces));
}

=head1 get_dna

Get the dna sequence from a genome/contig/beginning and end
my $dna = $minseed->get_dna($genome,$contig,$beg,$end);

=cut

sub get_dna {
	my ($self, $genome,$contig,$beg,$end) = @_;
	my $ln = ($end - $beg) + 1;
	return lc(substr($self->{$genome}->{dna_sequences}->{$contig}, ($beg-1),$ln));
}

=head1 between

Given an x, y, and z values are they inbetwen each other

=cut

sub between {
	my($self,$x,$y,$z) = @_;
	if ($x < $z) {
		return (($x <= $y) && ($y <= $z));
	} else {
		return (($x >= $y) && ($y >= $z));
	}
}


=head1 contigs_of

my @contigs=  $minseed->contigs_of($genome);

return a list of the contigs

=cut

sub contigs_of {	
	my ($self, $genome) = @_;
	if (!defined $self->{$genome}->{dna_sequences}) {$self->all_dna($genome)}
	return keys %{$self->{$genome}->{dna_sequences}};
}

=head1 contigs_lengths

Return a list of tuples of contigs and their lengths.

my @tples = $minseed->contigs_lengths($genome);
my ($contig, $length) = @$tple;

=cut

sub contigs_lengths {
	my ($self, $genome) = @_;
	if (!defined $self->{$genome}->{dna_sequences}) {$self->all_dna($genome)}
	my @ret;
	map {push @ret, [$_, length($self->{$genome}->{dna_sequences}->{$_})]} keys %{$self->{$genome}->{dna_sequences}};
	return @ret;
}

=head1 contigs_lengths



=head1 contig_ln 

my $len = $miseed->contig_ln($genome, $contig);

get the length of the contig

=cut


sub contig_ln {
	my ($self, $genome, $contig) = @_;
	if (!defined $self->{$genome}->{dna_sequences}) {$self->all_dna($genome)}
	if (!defined $self->{$genome}->{dna_sequences}->{$contig}) {die "$contig does not appear to be in $genome"}
	return length($self->{$genome}->{dna_sequences}->{$contig});
}



=head1 all_dna

Read all the dna sequences from a genome

=cut 

sub all_dna {
	my ($self, $genome) = @_;
	return undef if (!defined $self->{orgdir});
	$self->{$genome}->{dna_sequences}=$rob->read_fasta($self->{orgdir}."/$genome/contigs");
	return $self->{$genome}->{dna_sequences};
}

=head1 function_of

Print the function of a protein

 my $function = $minseed->function_of($peg);

NOTE: This is different than the typical function_of because there is not an array version

=cut

sub function_of {
	my ($self, $peg) = @_;
	return undef unless ($self->{orgdir});
	return $self->{function_of}->{$peg} if (defined $self->{function_of}->{$peg});
	if ($self->{$self->genome_of($peg)}->{is_deleted}->{$peg}) {print STDERR "WARNING: $peg is deleted\n"}
	my %functions;
	my $genome = $self->genome_of($peg);
	foreach my $annfile (qw[proposed_non_ff_functions proposed_functions assigned_functions]) {
		next unless (-e  $self->{orgdir}."/$genome/$annfile");
		open(IN, $self->{orgdir}."/$genome/$annfile") || die "can't open " . $self->{orgdir}."/$genome/$annfile";
		while (<IN>) {
			chomp;
			my @a=split /\t/;
			$functions{$a[0]}=$a[1];
		}
	}
	foreach my $allpeg ($self->pegs_of($self->genome_of($peg))) {
		$self->{function_of}->{$allpeg} = (defined  $functions{$allpeg}) ? $functions{$allpeg} : "hypothetical protein";
	}
	return $self->{function_of}->{$peg};
}

=head1 assign_function

Assign a function to a peg. This will set the function for that peg. If you provide a who or a why, this will also append to the annotations file. You should probably do this, but it is not essential for MinSeed does not use that file.

my $worked = $fig->assign_function($peg, $who, $function, $why);

=cut

sub assign_function {
	my ($self, $peg, $who, $function, $why) = @_;
	return 0 unless ($peg && $function);
	my $genome = $self->genome_of($peg);
	open(OUT,  ">>".$self->{orgdir}."/$genome/assigned_functions") || die "can't append to assigned_functions";
	print OUT "$peg\t$function\n";
	close OUT;

	if ($who || $why) {
		open(OUT, ">>".$self->{orgdir}."/$genome/annotations") || die "can't open to annotations";
		print OUT join("\n", $peg, time, $who, "Set function to", $function, $why, "//"), "\n";
		close OUT;
	}
	return 1;
}




=head1 get_translation

Get the translation of a protein

my $trans = $fig->get_translation($peg);

=cut

sub get_translation {
	my ($self, $peg) = @_;
	return $self->{'translation'}->{$peg} if ($self->{'translation'}->{$peg});
	if ($self->{'orgdir'}) {
		$self->{'translation'} = $rob->read_fasta($self->{'orgdir'}."/".$self->genome_of($peg)."/Features/peg/fasta");
		return $self->{'translation'}->{$peg};
	}
	else {
		return "";
	}
}

=head1 set_translation

Set the translation of a protein

my $worked = $fig->set_translation($peg, $translation);

=cut

sub set_translation {
	my ($self, $peg, $translation) = @_;
	$self->{'translation'}->{$peg}=$translation;
	my $fastaf = $self->{'orgdir'}."/".$self->genome_of($peg)."/Features/".$self->ftype($peg)."/fasta";
	unless (-e $fastaf) {die "Can not find the existing fasta file we think is at $fastaf"}
	open(TRA, ">>$fastaf") || die "can't append to $fastaf";
	print TRA ">$peg\n$translation\n";
	close TRA;
	return 1;
}
		
=head1 pegs_in_order

Get the genes in order along the genome. This is code stolen from elsewhere :)

my @pegs = $fig->pegs_in_order($genome);

=cut

sub pegs_in_order {
        my ($self, $genome) = @_;
        my @pegs = map  { $_->[0] }
        sort { ($a->[1] cmp $b->[1]) or ($a->[2] <=> $b->[2]) }
        map  { my $peg = $_;
                if (my $loc = $self->feature_location($peg) )
                {
                        my ($contig,$beg,$end) = $self->boundaries_of($loc);
                        [$peg,$contig,$self->min($beg,$end)];
                }
                else
                {
                        ();
                }
        }
        $self->pegs_of($genome);
        return @pegs;
}	


=head1 boundaries_of

get the boundaries of a location

my ($contig,$beg,$end) = $fig->boundaries_of($loc);

=cut

sub boundaries_of {
	my $self=shift;
	my($location) = (@_ == 1) ? $_[0] : $_[1];
	my($contigQ);

	if (defined($location))
	{
		my @exons = split(/\s*,\s*/,$location);
		my($contig, $beg, $end);

		if (($exons[0] =~ /^(\S+)_(\d+)_\d+$/))
		{
			($contig, $beg) = ($1,$2);
			$contigQ = quotemeta $contig;

			if ($exons[$#exons] =~ /^$contigQ\_\d+_(\d+)$/)
			{
				$end = $1;
				if ($beg > 0 && $end > 0)
				{
					my $strand = (($beg < $end) ? qq(+) : qq(-));
					return ($contig, $beg, $end, $strand);
				}
			}
		}
		Cluck("Could not parse loc=$location.") if T(0);
	}
	return ();
}

=head1 max

Find the maximum of two values

$max = $fig->max($a, $b);

=cut

sub max {
	my ($self, $x, $y) = @_;
	return ($x > $y) ? $x : $y;
}

=head1 min

Find the minimum of two values

$min = $fig->min($a, $b);

=cut

sub min {
	my ($self, $x, $y) = @_;
	return ($x < $y) ? $x : $y;
}

=head1 sims

Get the sims for a protein

my @sims = $fig->sims($peg, $max, $maxP, $select)

peg is the protein
max is the maximum number of pegs to return (default 10000)
maxP is the maximum E value (default 1e-5)
select:
	if undef or all then it will return all sims
	if figx will only return sims with a fig id

=cut



sub sims {
	my ($self, $peg, $max, $maxP, $select) = @_;
	$max     = $max ? $max : 10000;
	$maxP    = $maxP ? $maxP : 1.0e-5;

	if ($self->{'sims_genome'} != $self->genome_of($peg)) { 
		$self->load_sims_index($self->genome_of($peg)) ;
	}
	
	my $location = $self->{"sims_index"}->{$peg};
	return if (!$location || $location !~ /^(\d+),(\d+)$/); # this is not a valid location or the peg has no sims
	my($seek, $len) = ($1, $2); # where to read the file
	my $sims_txt = $self->read_block($self->{'sims_fh'}, $seek, $len); # read the file
	
	my @sims = grep { $_->[10] <= $maxP } map { chomp; [split(/\t/)] } @$sims_txt;

	if ($select eq "figx") { @sims = grep { index( $_->[1], 'fig') > -1 } @sims }

	return splice(@sims, 0, $max);

}

sub load_sims_index {
	my ($self, $genome) = @_;

	my $sims_file = $self->{orgdir}. "/$genome/expanded_similarities";
	my $index_file = "$sims_file.index";
	my $hash = {};

	my $tied = tie %$hash, 'DB_File', $index_file, O_RDONLY, 0666, $DB_BTREE;

#
# Set these even if failed so we don't keep trying to open and failing.
#
	$self->{"sims_index"} = $hash;
	$self->{"sims_tie"} = $tied;

	if (not $tied)
	{
		warn "Cannot tie sims index $index_file: $!\n";
	}

#
# open the sims file as well.
#

	$self->{'sims_fh'} = new FileHandle("<$sims_file");

	if (!$self->{"sims_fh"})
	{
		warn "Cannot open sims file $sims_file: $!\n";
	}
	$self->{'sims_genome'}=$genome;
}


# read a block from a file
sub read_block {
	my($self, $fh, $seek, $ln) = @_;
	my($piece,$readN);

	seek($fh,$seek,0);
	my @lines = ();

	$readN = read($fh,$piece,$ln);
	($readN == $ln) || print STDERR "could not read the block of sims at $seek for $ln characters; $readN actually read";
	return [ split( /\n/, $piece ) ];
}




=head1 same_strand

Given a bunch of locations in format contig_beg_end are they all on the same strand

my $flag = $fig->smae_strand($loc1, $loc2, ...);

=cut

sub same_strand {
	my ($self, @locs)=@_;
	my $rev; my $fwd;
	foreach my $loc (@locs) {
		my ($contig, $beg, $end, $strand) = $self->boundaries_of($loc);
		($strand eq "+") ? ($fwd = 1) : ($strand eq "-") ? ($rev = 1) : print STDERR "Not sure what $strand strand is \n";
	}
	($rev && !$fwd) ? return 1 :
		(!$rev && $fwd) ? return 1 :
		return 0;
}



=head1 next_available_protein_id

Get the next available protein id for this genome. When we delete pegs, we mark them as deleted, but don't delete them. We never re use an id, so this just gets the next available ID

my $id = $fig->next_available_id($genome);

=cut

sub next_available_id {
	my ($self, $genome) = @_;
	unless ($self->{'last_id'}->{$genome}) {
		$self->pegs_of($genome);
	}
	$self->{'last_id'}->{$genome}++;
	return $self->{'last_id'}->{$genome};
}


=head1 subsystems_for_peg

Return the list of subsystems and roles that this peg appears in.
Returns an array. Each item in the array is
a reference to a tuple of subsystem and role.

=cut

sub subsystems_for_peg {
	my ($self, $peg) = @_;
	my $genome = $self->genome_of($peg);
	unless ($self->{'subsystems'}->{$genome}) {
		my $ssfile;
		if (-e $self->{orgdir} . "/$genome/Subsystems/bindings_updated") {
			$ssfile = $self->{orgdir} . "/$genome/Subsystems/bindings_updated";
		}
		elsif (-e $self->{orgdir} . "/$genome/Subsystems/bindings") {
			$ssfile = $self->{orgdir} . "/$genome/Subsystems/bindings";
		}
		else {
			print STDERR "No subsystem bindings found for $genome\n";
			$self->{'subsystems'}->{$genome}={};
			return [];
		}
		open(IN, $ssfile) || die "Can't open $ssfile";
		while (<IN>) {
			chomp;
			my ($ss, $role, $peg)=split /\t/;
			push @{ $self->{'subsystems'}->{$genome}->{$peg} }, [$ss, $role];
		}
		close IN;
	}

	if (defined $self->{'subsystems'}->{$genome}->{$peg}) {
		return $self->{'subsystems'}->{$genome}->{$peg};
	} else {
		return [];
	}
}

=head1 ss_class

Get the two level classification for the subsystems
my $class = $data->ss_class($subsys)

Returns a reference to an array or ["", ""]

=cut

sub ss_class {
	my ($self, $ss)=@_;
	if ($self->{classification}->{$ss}) {
		return $self->{classification}->{$ss};
	}

	my @genomes = $self->genomes();
	open(IN, $self->{orgdir} . "/".$genomes[0]. "/Subsystems/classification") || die "No classification. Did you run /data/Vibrio/Rob/Reannotate_subsystems/add_ss_class.pl";
	while (<IN>) {
		chomp;
		my @a=split /\t/;
		$self->{classification}->{$a[2]}=[$a[0], $a[1]];
	}
	close IN;
	if (!$self->{classification}->{$ss}) {
		$self->{classification}->{$ss}=["", ""];
	}
	return $self->{classification}->{$ss};
}



=head1 is_hypothetical

Returns true if the function is hypothetical

if (!$data->is_hypothetical("a text string for the function")) { ... }

=cut

sub is_hypothetical {
	my ($self, $x)=@_;
	if (! $x)                             { return 1 }
	if ($x =~ /lmo\d+ protein/i)          { return 1 }
	if ($x =~ /hypoth/i)                  { return 1 }
	if ($x =~ /conserved protein/i)       { return 1 }
	if ($x =~ /gene product/i)            { return 1 }
	if ($x =~ /interpro/i)                { return 1 }
	if ($x =~ /B[sl][lr]\d/i)             { return 1 }
	if ($x =~ /^U\d/)                     { return 1 }
	if ($x =~ /^orf[^_]/i)                { return 1 }
	if ($x =~ /uncharacterized/i)         { return 1 }
	if ($x =~ /pseudogene/i)              { return 1 }
	if ($x =~ /^predicted/i)              { return 1 }
	if ($x =~ /AGR_/)                     { return 1 }
	if ($x =~ /similar to/i)              { return 1 }
	if ($x =~ /similarity/i)              { return 1 }
	if ($x =~ /glimmer/i)                 { return 1 }
	if ($x =~ /unknown/i)                 { return 1 }
	if (($x =~ /domain/i) ||
			($x =~ /^y[a-z]{2,4}\b/i) ||
			($x =~ /complete/i) ||
			($x =~ /ensang/i) ||
			($x =~ /unnamed/i) ||
			($x =~ /EG:/) ||
			($x =~ /orf\d+/i) ||
			($x =~ /RIKEN/) ||
			($x =~ /Expressed/i) ||
			($x =~ /[a-zA-Z]{2,3}\|/) ||
			($x =~ /predicted by Psort/) ||
			($x =~ /^bh\d+/i) ||
			($x =~ /cds_/i) ||
			($x =~ /^[a-z]{2,3}\d+[^:\+\-0-9]/i) ||
			($x =~ /similar to/i) ||
			($x =~ / identi/i) ||
			($x =~ /ortholog of/i) ||
			($x eq "Phage protein") ||
			($x =~ /structural feature/i))    { return 1 }
	return 0;
}





1;
