#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#

#
# This is a SAS component.
#

package SAP;

    use strict;
    use ERDB;
    use Tracer;
    use SeedUtils;
    use ServerThing;

=head1 Sapling Server Function Object

This file contains the functions and utilities used by the Sapling Server
(B<sap_server.cgi>). The various methods listed in the sections below represent
function calls direct to the server. These all have a signature similar to the
following.

    my $results = $sapObject->function_name($args);

where C<$sapObject> is an object created by this module,
C<$args> is a parameter structure, and C<function_name> is the Sapling
Server function name. The output $results is a scalar, generally a hash
reference, but sometimes a string or a list reference.

=head2 Location Strings

Several methods deal with gene locations. Location information from the Sapling
server is expressed as I<location strings>. A location string consists of a
contig ID (which includes the genome ID), an underscore, a starting location, a
strand indicator (C<+> or C<->), and a length. The first location on the contig
is C<1>.

For example, C<100226.1:NC_003888_3766170+612> indicates contig C<NC_003888> in
genome C<100226.1> (I<Streptomyces coelicolor A3(2)>) beginning at location
3766170 and proceeding forward on the plus strand for 612 bases.

=head2 Constructor

Use

    my $sapObject = SAPserver->new();

to create a new sapling server function object. The server function object
is used to invoke the L</Primary Methods> listed below. See L<SAPserver> for
more information on how to create this object and the options available.

=cut

#
# Actually, if you are using SAP.pm, you should do SAP->new(), not SAPserver->new()
# That comment above is for the benefit of the pod doc stuff on how to use SAPserver
# that is generated from this file.
#

sub new {
    my ($class, $sap) = @_;
    # Create the sapling object.
    if (! defined $sap) {
        $sap = ERDB::GetDatabase('Sapling');
    }
    # Create the server object.
    my $retVal = { db => $sap };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

sub _set_memcache
{
    my($self, $mcache) = @_;
    $self->{memcache} = $mcache;
}



=head1 Primary Methods

=head2 Server Utility Methods

You will not use the methods in this section very often. Some are used by the
server framework for maintenance and control purposes (L</methods>), while others
(L</query> and L</get>) provide access to data in the database in case you need
data not available from one of the standard methods.

=head3 methods

    my $methodList =        $sapObject->methods();

Return a reference to a list of the methods allowed on this object.

=cut

use constant METHODS => [qw(
                            all_complexes
                            all_experiments
                            all_features
                            all_figfams
                            all_genomes
                            all_models
                            all_proteins
                            all_reactions
                            all_roles_used_in_models
                            all_subsystems
                            atomic_regulon_vectors
                            atomic_regulons
                            classification_of
                            close_genomes
                            clusters_containing
                            co_occurrence_evidence
                            compared_regions
                            complex_data
                            conserved_in_neighborhood
                            contig_lengths
                            contig_sequences
                            coregulated_correspondence
                            coregulated_fids
                            coupled_reactions
                            discriminating_figfams
                            dlits_for_ids
                            equiv_ids_for_sequences
                            equiv_precise_assertions
                            equiv_sequence_assertions
                            equiv_sequence_ids
                            exists
                            experiment_fid_levels
                            experiment_regulon_levels
                            expressed_genomes
                            feature_assignments
                            fid_correspondences
                            fid_experiments
                            fid_locations
                            fid_map_for_genome
                            fid_possibly_truncated
                            fid_vectors
                            fids_expressed_in_range
                            fids_to_ids
                            fids_to_proteins
                            fids_to_regulons
                            fids_with_evidence_code
                            fids_with_evidence_codes
                            figfam_fids
                            figfam_fids_batch
                            figfam_function
                            find_closest_genes
                            genes_in_region
                            gene_correspondence_map
                            genome_contigs
                            genome_contig_md5s
                            genome_data
                            genome_domain
                            genome_fid_md5s
                            genome_experiments
                            genome_experiment_levels
                            genome_figfams
                            genome_ids
                            genome_metrics
                            genome_names
                            genomes_to_subsystems
                            genomes_by_md5
                            get
                            get_subsystems
                            ids_in_subsystems
                            ids_to_annotations
                            ids_to_assertions
                            ids_to_data
                            ids_to_fids
                            ids_to_figfams
                            ids_to_functions
                            ids_to_genomes
                            ids_to_lengths
                            ids_to_publications
                            ids_to_sequences
                            ids_to_subsystems
                            intergenic_regions
                            is_in_subsystem
                            is_in_subsystem_with
                            is_prokaryotic
                            locs_to_dna
                            make_runs
                            mapped_genomes
                            models_to_reactions
                            occ_of_role
                            otu_members
                            pairsets
                            pegs_implementing_roles
                            pegs_in_subsystem
                            pegs_in_subsystems
                            pegs_in_variants
                            proteins_to_fids
                            query
                            reaction_neighbors
                            reaction_path
                            reaction_strings
                            reactions_to_complexes
                            reactions_to_roles
                            regulons_to_fids
                            related_clusters
                            related_figfams
                            representative
                            representative_genomes
                            role_neighbors
                            role_reactions
                            roles_exist_in_subsystem
                            roles_to_complexes
                            roles_to_figfams
                            roles_to_proteins
                            roles_to_subsystems
                            rows_of_subsystems
                            scenario_names
                            select
                            submit_gene_correspondence
                            subsystem_data
                            subsystem_genomes
                            subsystem_names
                            subsystem_roles
                            subsystem_spreadsheet
                            subsystem_type
                            subsystems_for_role
                            taxonomy_of
                            upstream
                        )];

sub methods {
    # Get the parameters.
    my ($self) = @_;
    # Return the result.
    return METHODS;
}


=head3 exists

    my $idHash =            $sapObject->exists({
                                -type => 'Genome',
                                -ids => [$id1, $id2, ...]
                            });

Return a hash indicating which of the specified objects of the given type exist
in the database. This method is used as a general mechanism for finding what
exists and what doesn't exist when you know the ID. In particular, you can use
it to check for the presence or absence of subsystems, genomes, features,
or FIGfams.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -type

The type of object whose existence is being queried. The type specification is
case-insensitive: C<genome> and C<Genome> are treated the same. The permissible
types are

=over 12

=item Genome

Genomes, identified by taxon ID: C<100226.1>, C<83333.1>, C<360108.3>

=item Feature

Features (genes), identified by FIG ID: C<fig|100226.1.peg.3361>, C<fig|360108.3.rna.4>

=item Subsystem

Subsystem, identified by subsystem name: C<Arginine biosynthesis extended>

=item FIGfam

FIGfam protein family, identified by ID: C<FIG000171>, C<FIG001501>

=back

=item -ids

Reference to a list of identifiers for objects of the specified type.

=back

=item RETURN

Returns a reference to a hash keyed by ID. For each incoming ID, it maps
to C<1> if an object of the specified type with that ID exists, else C<0>.

    $idHash = { $id1 => $flag1, $id2 => $flag2, ... };

=back

=cut

use constant EXIST_OBJECT_TYPES => { genome => 'Genome', feature => 'Feature',
                                     subsystem => 'Subsystem', figfam => 'Family'
                                   };

sub exists {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of identifiers.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Get the object type.
    my $type = $args->{-type};
    Confess("No -type parameter specified.") if ! defined $type;
    my $objectType = EXIST_OBJECT_TYPES->{lc $type};
    Confess("Invalid object type \"$type\".") if ! defined $objectType;
    # Declare the return variable.
    my $retVal = {};
    # Loop through the identifiers, checking existence.
    for my $id (@$ids) {
        $retVal->{$id} = ($sap->Exists($objectType, $id) ? 1 : 0);
    }
    # Return the result.
    return $retVal;
}


=head3 get

    my $hashList =          $sapObject->get({
                                -objects => $objectNameString,
                                -filter => { $label1 => $criterion1, $label2 => $criterion2, ... },
                                -limit => $maxRows,
                                -fields => { $label1 => $name1, $label2 => $name2, ... },
                                -multiples => 'list',
                                -firstOnly => 1
                            });

Query the Sapling database. This is a variant of the L</query> method in
which a certain amount of power is sacrificed for ease of use. Instead of
a full-blown filter clause, the caller specifies a filter hash that maps
field identifiers to values.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -objects

The object name string listing all the entities and relationships in the
query. See L<ERDB/Object Name List> for more details.

=item -filter (optional)

Reference to a hash that maps field identifiers in L<ERDB/Standard Field Name Format>
to criteria. A criterion is either an object or scalar value (which is asserted
as the value of the field), a 2-tuple consisting of a relational operator and
a value (which is asserted to be in the appropriate relation to the field), or a
sub-list consisting of the word C<IN> and two or more values (which asserts that
the field has one of the listed values). A record satisfies the filter if it satisfies
all the criteria in the hash.

=item -limit (optional)

Maximum number of rows to return for this query. The default is no limit.

=item -fields (optional)

Reference to a hash mapping field identifiers to field names. In this case,
the field identifier is a field name in L<ERDB/Standard Field Name Format>
and the field name is the key value that will be used for the field in the
returned result hashes. If this parameter is omitted, then instead of a
returning the results, this method will return a count of the number of
records found.

=item -multiples (optional)

Rule for handling field values in the result hashes. The default option is
C<smart>, which maps single-valued fields to scalars and multi-valued fields
to list references. If C<primary> is specified, then all fields are mapped
to scalars-- only the first value of a multi-valued field is retained. If
C<list> is specified, then all fields are mapped to lists.

=item -firstOnly (optional)

If TRUE, only the first result will be returned. In this case, the return
value will be a hash reference instead of a list of hash references. The
default is FALSE.

=back

=item RETURN

Returns a reference to a list of hashes. Each hash represents a single record in
the result set, and maps the output field names to the field values for that
record. Note that if a field is multi-valued, it will be represented as a
list reference.

    $hashList = [{ $label1 => $row1value1, $label2 => $row1value2, ... },
                 { $label1 => $row2value1, $label2 => $row2value2, ... },
                 ... ];

=back

=cut

sub get {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the filter hash, the flags, and the limit. All of these
    # are optional.
    my $filter = $args->{-filter} || {};
    my $limit = $args->{-limit} || 0;
    my $multiples = $args->{-multiples} || 'smart';
    my $firstOnly = $args->{-firstOnly} || 0;
    # Initialize the return variable. In first-only mode, it starts
    # undefined; otherwise, it's an empty list.
    my $retVal = ($firstOnly ? undef : []);
    # Get the object name string and the result field hash.
    my $objects = $args->{-objects};
    my $fields = $args->{-fields};
    # Insure the object name list is present.
    if (! $objects) {
        Confess("Object name string not specified.");
    } else {
        # Get the default object name from the object name list.
        my ($defaultObject) = split m/\s+/, $objects, 2;
        # We'll build the filter elements and the parameter list in
        # these lists. The filter string will be formed by ANDing
        # together the filter elements.
        my (@filters, @parms);
        # Loop through the filter hash.
        for my $filterField (keys %$filter) {
            # Compute the field type.
            my $fieldType = $sap->FieldType($filterField);
            # Get this field's criterion.
            my $criterion = $filter->{$filterField};
            # If the criterion is a not a list, make it one.
            if (! defined $criterion) {
                Confess("Invalid (missing) criterion for field \"$filterField\".");
            } elsif (ref $criterion ne 'ARRAY') {
                $criterion = ['=', $criterion];
            }
            # Determine the criterion type.
            if ($criterion->[0] eq 'LIKE') {
                # For a LIKE, we don't convert the value.
                push @filters, "$filterField LIKE ?";
                push @parms, $criterion->[1];
            } elsif ($criterion->[0] eq 'IN') {
                # For an IN, we have to deal with multiple field values. We'll
                # stash one question mark per value in this list.
                my @marks;
                # Process the criterion elements.
                for (my $i = 1; $i < scalar @$criterion; $i++) {
                    push @marks, "?";
                    push @parms, $fieldType->encode($criterion->[$i]);
                }
                # Form the filter element from the collected marks.
                push @filters, "$filterField IN (" . join(", ", @marks) . ")";
            } else {
                # here we have a normal comparison.
                push @filters, "$filterField $criterion->[0] ?";
                push @parms, $fieldType->encode($criterion->[1]);
            }
        }
        # Create the filter string.
        my $filterString = join(" AND ", @filters);
        # Add the limit clause.
        if ($firstOnly) {
            $filterString .= " LIMIT 1";
        } elsif ($limit > 0) {
            $filterString .= " LIMIT $limit";
        }
        # Is this a query or a count?
        if (! defined $fields) {
            # It's a count. Do a GetCount call.
            $retVal = $sap->GetCount($objects, $filterString, \@parms);
        } else {
            # Here we have a real query. Now we run it.
            my $query = $sap->Get($objects, $filterString, \@parms);
            # Loop through the results.
            while (my $record = $query->Fetch()) {
                # Create the result hash for this record.
                my %results;
                # Loop through the fields.
                for my $outputField (keys %$fields) {
                    # Get the value.
                    my @values = $record->Value($outputField);
                    # Get the output field name.
                    my $outputName = $fields->{$outputField};
                    # Process according to the output type.
                    if ($multiples eq 'list') {
                        $results{$outputName} = \@values;
                    } elsif ($multiples eq 'primary' || scalar(@values) == 1) {
                        $results{$outputName} = $values[0];
                    } else {
                        $results{$outputName} = \@values;
                    }
                }
                # Add the result hash to the output. In first-only mode,
                # we store it; otherwise, we push it in.
                if ($firstOnly) {
                    $retVal = \%results;
                } else {
                    push @$retVal, \%results;
                }
            }
        }
    }
    # Return the result.
    return $retVal;
}


=head3 query

    my $rowList =           $sapObject->query({
                                -objects => $objectNameString,
                                -filterString => $whereString,
                                -limit => $maxRows,
                                -parameters => [$parm1, $parm2, ...],
                                -fields => [$name1, $name2, ...]
                            });

This method queries the Sapling database and returns a reference to a list of
lists. The query is specified in the form of an object name string, a filter
string, an optional list of parameter values, and a list of desired output
fields. The result document can be thought of as a two-dimensional array, with
each row being a record returned by the query and each column representing an
output field.

This function buys a great deal of flexibility as the cost of ease of use.
Before attempting to formulate a query, you will need to look at the
L<ERDB> documentation.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -objects

The object name string listing all the entities and relationships in the
query. See L<ERDB/Object Name List> for more details.

=item -filterString

The filter string for the query. It cannot contain a C<LIMIT> clause, but
can otherwise be anything described in L<ERDB/Filter Clause>.

=item -limit (optional)

Maximum number of rows to return for this query. The default is C<1000>. To
make an unlimited query, specify C<none>.

=item -parameters (optional)

Reference to a list of parameter values. These should be numbers or strings,
and are substituted for any parameter marks in the query on a one-for-one
basis. See also L<ERDB/Parameter List>.

=item -fields

Reference to a list containing the names of the desired output fields.

=back

=item RETURN

Returns a reference to a list of lists. Each row corresponds to a database
result row, and each column corresponds to one of the incoming output fields.
Note that some fields contain complex PERL data structures, and fields that
are multi-valued will contain sub-lists.

    $rowList = [[$row1field1, $row1field2, ...],
                [$row2field1, $row2field2, ...],
                [$row3field1, $row3field2, ...],
                ... ];

=back

=cut

sub query {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the query parameters.
    my $objects = $args->{-objects} || '';
    my $parms = $args->{-parameters} || [];
    my $filter = $args->{-filterString} || '';
    my $limitNumber = $args->{-limit} || 1000;
    my $fields = $args->{-fields} || [];
    # If this is an unlimited query, set the limit number to 0.
    if ($limitNumber eq 'none') {
        $limitNumber = 0;
    }
    # Declare the return variable.
    my @retVal;
    # Load the query console so we can get its help.
    require ERDBQueryConsole;
    # Create the console object. The user is allowed unlimited queries
    # because we encourage limits and we're hoping no one goes crazy.
    # The return data is meant to be raw rather than HTML.
    Trace("Submitting query.") if T(3);
    my $console = ERDBQueryConsole->new($self->{db}, secure => 1, raw => 1);
    # Try to submit the query.
    my $ok = $console->Submit($objects, $filter, $parms, $fields, $limitNumber);
    # Only proceed if there's no error.
    if (! $ok) {
        die $console->Messages();
    } else {
        Trace("Processing query results.") if T(3);
        # Loop through the result rows.
        while (my @row = $console->GetRow()) {
            push @retVal, \@row;
        }
    }
    # Return the result.
    return \@retVal;
}

=head3 select

    my $listList =          $sapObject->select({
                                -path => $objectNameString,
                                -filter => { $field1 => $list1, $field2 => $list2, ... },
                                -fields => [$fieldA, $fieldB, ... ],
                                -limit => $maxRows,
                                -multiples => 'list'
                            });

Query the Sapling database. This is a variant of the L</get> method in
which a further amount of power is sacrificed for ease of use. The
return is a list of lists, and the criteria are always in the form of
lists of possible values.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -path

The object name string listing all the entities and relationships in the
query. See L<ERDB/Object Name List> for more details.

=item -filter (optional)

Reference to a hash that maps field identifiers in L<ERDB/Standard Field Name Format>
to lists of permissible values. A record matches the filter if the field value
matches at least one element of the list.

=item -fields

Reference to a list of field names in L<ERDB/Standard Field Name Format>.

=item -limit (optional)

Maximum number of rows to return for this query. The default is no limit.

=item -multiples (optional)

Rule for handling field values in the result hashes. The default option is
C<smart>, which maps single-valued fields to scalars and multi-valued fields
to list references. If C<primary> is specified, then all fields are mapped
to scalars-- only the first value of a multi-valued field is retained. If
C<list> is specified, then all fields are mapped to lists.

=back

=item RETURN

Returns a reference to a list of lists. Each sub-list represents a single record
in the result set, and contains the field values in the order the fields were
lists in the C<-fields> parameter. Note that if a field is multi-valued, it will
be represented as a list reference.

    $listList = [[$row1value1, $row1value2, ... ], [$row2value1, $row2value2, ...], ... ];

=back

=cut

sub select {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the filter hash, the flags, and the limit. All of these
    # are optional.
    my $filter = $args->{-filter} || {};
    my $limit = $args->{-limit} || 0;
    my $multiples = $args->{-multiples} || 'smart';
    # Initialize the return variable.
    my $retVal = [];
    # Get the object name string and the result field list.
    my $objects = $args->{-path};
    my $fields = ServerThing::GetIdList(-fields => $args);
    # Insure the object name list is present.
    if (! $objects) {
        Confess("Object name string not specified.");
    } else {
        # Get the default object name from the object name list.
        my ($defaultObject) = split m/\s+/, $objects, 2;
        # We'll build the filter elements and the parameter list in
        # these lists. The filter string will be formed by ANDing
        # together the filter elements.
        my (@filters, @parms);
        # Loop through the filter hash.
        for my $filterField (keys %$filter) {
            # Compute the field type.
            my $fieldType = $sap->FieldType($filterField);
            # Get this field's criterion.
            my $criterion = $filter->{$filterField};
            # Insure the criterion exists.
            if (! defined $criterion) {
                Confess("Invalid (missing) criterion for field \"$filterField\".");
            } elsif (ref $criterion ne 'ARRAY') {
                # Here we have a scalar criterion. It is encoded as
                # an equality clause.
                push @parms, $fieldType->encode($criterion);
                push @filters, "$filterField = ?";
            } else {
                # Here we have to deal with multiple field values. We'll
                # stash one question mark per value in this list.
                my @marks;
                # Process the criterion elements.
                for (my $i = 0; $i < scalar @$criterion; $i++) {
                    push @marks, "?";
                    push @parms, $fieldType->encode($criterion->[$i]);
                }
                # Form the filter element from the collected marks.
                push @filters, "$filterField IN (" . join(", ", @marks) . ")";
            }
        }
        # Create the filter string.
        my $filterString = join(" AND ", @filters);
        # Add the limit clause.
        if ($limit > 0) {
            $filterString .= " LIMIT $limit";
        }
        # Run the query.
        my $query = $sap->Get($objects, $filterString, \@parms);
        # Loop through the results.
        while (my $record = $query->Fetch()) {
            # Create the result list for this record.
            my @results;
            # Loop through the fields.
            for my $outputField (@$fields) {
                # Get the value.
                my @values = $record->Value($outputField);
                # Process according to the output type.
                if ($multiples eq 'list') {
                    push @results, \@values;
                } elsif ($multiples eq 'primary' || scalar(@values) == 1) {
                    push @results, $values[0];
                } else {
                    push @results, \@values;
                }
            }
            # Add the result list to the output. In first-only mode,
            # we store it; otherwise, we push it in.
            push @$retVal, \@results;
        }
    }
    # Return the result.
    return $retVal;
}

=head2 Annotation and Assertion Data Methods

=head3 equiv_precise_assertions

    my $idHash =            $sapObject->equiv_precise_assertions({
                                -ids => [$id1, $id2, ...]
                            });

Return the assertions for all genes in the database that match the
identified gene. The gene can be specified by any prefixed gene
identifier (e.g. C<uni|AYQ44>, C<gi|85841784>, or
C<fig|360108.3.peg.1041>).

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of gene identifiers.

=back

For backward compatibility, the parameter can also be a reference to a list
of gene identifiers.

=item RETURN

Returns a reference to a hash that maps each incoming ID to a list of 4-tuples.
Each 4-tuple contains (0) an identifier that is for the same gene as the input
identifier, (1) the asserted function of that identifier, (2) the source of
the assertion, and (3) a flag that is TRUE if the assertion is by a human expert.

    $idHash = { $id1 => [$otherID1, $function1, $source1, $flag1],
                $id2 => [$otherID2, $function2, $source2, $flag2],
                ... };

In backward-compatibility mode, returns a reference to a list of 2-tuples. Each
2-tuple consists of an incoming ID and the list of 4-tuples with the asserted
function information.

=back

=cut

sub equiv_precise_assertions {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Check for backward compatibility mode.
    my $backwardMode = 0;
    if (ref $args ne 'HASH') {
        $args = { -ids => $args };
        $backwardMode = 1;
    }
    # Get the list of IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    foreach my $id (@$ids) {
        my @resultRows = $sap->GetAll("Identifier HasAssertionFrom Source",
                                      'Identifier(id) = ? ',
                                      [$id], [qw(Identifier(id)
                                                 HasAssertionFrom(function)
                                                 Source(id)
                                                 HasAssertionFrom(expert))]);
        $retVal->{$id} = \@resultRows;
    }
    # Check for backward-compatibility mode.
    if ($backwardMode) {
        # Convert the hash to a list of 2-tuples.
        my @outList = map { [$_, $retVal->{$_}] } @$ids;
        $retVal = \@outList;
    }
    # Return the result.
    return $retVal;
}

=head3 equiv_sequence_assertions

    my $idHash =            $sapObject->equiv_sequence_assertions({
                                -ids => [$id1, $id2, ...]
                            });

Return the assertions for all genes in the database that match the
identified protein sequences. A protein sequence can be identified by a
protein MD5 code or any prefixed gene identifier (e.g. C<uni|AYQ44>,
C<gi|85841784>, or C<fig|360108.3.peg.1041>).

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of protein identifiers. Each identifier should be a prefixed
gene identifier or the (optionally) prefixed MD5 of a protein sequence.

=back

=item RETURN

Returns a reference to a hash mapping each incoming protein identifier to a list
of 5-tuples, consisting of (0) an identifier that is sequence-equivalent to the
input identifier, (1) the asserted function of that identifier, (2) the source
of the assertion, (3) a flag that is TRUE if the assertion is by an expert, and
(4) the name of the genome relevant to the identifer (if any).

    $idHash = { $id1 => [$otherID1, $function1, $source1, $flag1],
                $id2 => [$otherID2, $function2, $source2, $flag2],
                ... };

=back

=cut

sub equiv_sequence_assertions {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Convert a list to a hash.
    if (ref $args ne 'HASH') {
        $args = { -ids => $args };
    }
    # Declare the return variable.
    my $retVal = {};
    # Get the list of IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs in the list.
    for my $id (@$ids) {
        # This hash will contain a list of the relevant protein sequence IDs.
        my %prots;
        # We'll put our assertions found in here.
        my @results;
        # Determine the ID type.
        if (my $prot = $sap->IsProteinID($id)) {
            # Here we have a protein sequence MD5 ID. In this case, we just
            # strip the prefix to get a Sapling protein sequence ID.
            $prots{$prot} = 1;
        } else {
            # Here we have a gene ID. Start by asking for all of the
            # protein sequences it identifies directly.
            my @prots = $sap->GetFlat("Identifier Names ProteinSequence",
                                      'Identifier(id) = ?', [$id],
                                      'ProteinSequence(id)');
            # Add the ones it identifies through a feature.
            push @prots, $sap->GetFlat("Identifier Identifies Feature Produces ProteinSequence",
                                       'Identifier(id) = ?', [$id],
                                       'ProteinSequence(id)');
            # Put all the proteins found in the hash.
            for my $prot (@prots) {
                $prots{$prot} = 1;
            }
        }
        # Loop through the protein sequences, finding assertions. For each
        # protein, we make two queries. Note that we expect the number of
        # protein sequences to be small, despite the large amount of work
        # performed above.
        for my $prot (sort keys %prots) {
            # Get the assertions on the protein's identifiers.
            @results = $sap->GetAll("ProteinSequence IsNamedBy Identifier HasAssertionFrom Source",
                                    "ProteinSequence(id) = ?", [$prot],
                                    [qw(Identifier(id) HasAssertionFrom(function)
                                        Source(id) HasAssertionFrom(expert))]);
            # Add the assertions on the identifiers for the protein's features.
            push @results, $sap->GetAll("ProteinSequence IsProteinFor Feature IsIdentifiedBy Identifier HasAssertionFrom Source AND Feature IsOwnedBy Genome",
                                        "ProteinSequence(id) = ?", [$prot],
                                        [qw(Identifier(id) HasAssertionFrom(function)
                                           Source(id) HasAssertionFrom(expert)
                                           Genome(scientific-name))]);
        }
        # If we found results, put them in the return object.
        Trace(scalar(@results) . " results found for $id.") if T(3);
        $retVal->{$id} = \@results;
    }
    # Return the result.
    return $retVal;
}

=head3 feature_assignments

    my $featureHash =       $sapObject->feature_assignments({
                                -genome => $genomeID,
                                -type => 'peg',
                                -hypothetical => 1
                            });

Return all features of the specified type for the specified genome along
with their assignments.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -genome

ID of the genome whose features are desired.

=item -type (optional)

If specified, the type of feature desired (C<peg>, C<rna>, etc.). If omitted,
all features will be returned.

=item -hypothetical (optional)

If C<1>, only hypothetical genes will be returned; if C<0>, only non-hypothetical
genes will be returned. If undefined or not specified, all genes will be
returned.

=back

=item RETURN

Returns a hash mapping the ID of each feature in the specified genome to
its assignment.

    $featureHash = { $fid1 => $function1, $fid2 => $function2, ... };

=back

=cut

sub feature_assignments {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the ID of the desired genome.
    my $genomeID = $args->{-genome};
    if (! $genomeID) {
        Confess("No genome ID specified for feature_assignments.");
    } else {
        # Start the feature ID search pattern.
        my $pattern = "fig|$genomeID.";
        # Add the feature type (if any).
        my $type = $args->{-type};
        if ($type) {
            $pattern .= "$type.";
        }
        # Append a wild card so we match everything in the genome.
        $pattern .= "%";
        # Ask for all the features.
        my @rows = $sap->GetAll("Feature", 'Feature(id) LIKE ?', [$pattern],
                                [qw(id function)]);
        # Are we checking for hypotheticals?
        my $hypothetical = $args->{-hypothetical};
        if (! defined $hypothetical) {
            # No. Put everything in the return hash.
            $retVal = { map { $_->[0] => $_->[1] } @rows };
        } else {
            # Yes. We need to check every functional assignment.
            for my $row (@rows) {
                # Get the ID and assignment.
                my ($fid, $assignment) = @$row;
                # Check to see if the assignment is hypothetical.
                my $rowHypothetical = (hypo($assignment) ? 1 : 0);
                Trace("Assignment \"$assignment\" has hypo = $rowHypothetical.") if T(3);
                # Include it if it matches the criterion specified by the caller.
                if ($rowHypothetical == $hypothetical) {
                    $retVal->{$fid} = $assignment;
                }
            }
        }
    }
    # Return the result.
    return $retVal;
}


=head3 ids_to_assertions

    my $idHash =            $sapObject->ids_to_assertions({
                                -ids => [$id1, $id2, ...]
                            });

Return the assertions associated with each prefixed ID.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of prefixed feature IDs (e.g. C<gi|17017961>,
C<NP_625335.1>, C<fig|360108.3.peg.1041>). The assertions associated with
each particular identifier will be returned. In this case, there will be
no processing for equivalent IDs. For that, you should use
L<equiv_sequence_assertions> or L<equiv_precise_assertions>.

=back

=item RETURN

Returns a reference to a hash mapping every incoming ID to a list of
3-tuples, each consisting of (0) an asserted function, (1)
the source of the assertion, and (2) a flag that is TRUE if the assertion
was made by an expert.

    $idHash = { $id1 => [[$assertion1a, $source1a, $expert1a],
                         [$assertion1b, $source1b, $expert1b], ...],
                $id2 => [[$assertion2a, $source2a, $expert2a],
                         [$assertion2b, $source2b, $expert2b], ...],
                ... };

=back

=cut

sub ids_to_assertions {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the feature ID list.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs.
    for my $id (@$ids) {
        # Get the assertion data.
        my @resultRows = $sap->GetAll("Identifier HasAssertionFrom Source",
                                      'Identifier(id) = ? ',
                                      [$id], [qw(HasAssertionFrom(function)
                                                 Source(id)
                                                 HasAssertionFrom(expert))]);
        # Store it in the return hash.
        $retVal->{$id} = \@resultRows;
    }
    # Return the results.
    return $retVal;
}


=head3 ids_to_annotations

    my $idHash =            $sapObject->ids_to_annotations({
                                -ids => [$id1, $id2, ...]
                            });

Return the annotations associated with each prefixed ID. Annotations are
comments attached to each feature (gene), and include past functional
assignments as well as more general information.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of feature IDs.

=item -source (optional)

Database source of the IDs specified-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
Use C<mixed> to allow mixed ID types (though this may cause problems when the same
ID has different meanings in different databases). Use C<prefixed> to allow IDs with
prefixing indicating the ID type (e.g. C<uni|P00934> for a UniProt ID, C<gi|135813> for
an NCBI identifier, and so forth). The default is C<SEED>.

=back

=item RETURN

Returns a reference to a hash mapping every incoming ID to a list of
3-tuples, each consisting of (0) annotation text, (1) the name of the
annotator, and (2) the timestamp of the annotation (as a number of seconds
since the epoch).

    $idHash = { $id1 => [[$annotation1a, $name1a, $time1a],
                         [$annotation1b, $name1b, $time1b], ...],
                $id2 => [[$annotation2a, $name2a, $time2a],
                         [$annotation2b, $name2b, $time2b], ...],
                ... };

=back

=cut

sub ids_to_annotations {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the feature ID list.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Get the desired ID type.
    my $source = $args->{-source};
    # Build the filter string, object name list, and parameter value list prefix.
    my ($objects, $filter, @parms) = $sap->ComputeFeatureFilter($source);
    # Loop through the incoming features.
    for my $id (@$ids) {
        # Get the annotation data.
        my @resultRows = $sap->GetAll("$objects IsAnnotatedBy Annotation",
                                      $filter, [@parms, $id],
                                      [qw(Annotation(comment)
                                          Annotation(annotator)
                                          Annotation(annotation-time))]);
        # Store it in the return hash.
        $retVal->{$id} = \@resultRows;
    }
    # Return the results.
    return $retVal;
}


=head3 ids_to_functions

    my $featureHash =       $sapObject->ids_to_functions({
                                -ids => [$id1, $id2, ...],
                                -source => 'CMR'
                                -genome => $genome
                            });

Return the functional assignment for each feature in the incoming list.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of feature IDs.

=item -source (optional)

Database source of the IDs specified-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
Use C<mixed> to allow mixed ID types (though this may cause problems when the same
ID has different meanings in different databases). Use C<prefixed> to allow IDs with
prefixing indicating the ID type (e.g. C<uni|P00934> for a UniProt ID, C<gi|135813> for
an NCBI identifier, and so forth). The default is C<SEED>.

=item -genome (optional)

ID of a specific genome. If specified, results will only be returned for genes in the
specified genome. The default is to return genes for all genomes.

=back

=item RETURN

Returns a reference to a hash mapping each feature ID to the feature's current
functional assignment. Features that do not exist in the database will not be
present in the hash. For IDs that correspond to multiple features, only one
functional assignment will be returned.

    $featureHash = { $id1 => $function1,
                     $id2 => $function2,
                     ...};

=back

=cut

sub ids_to_functions {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the feature ID list.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Get the desired ID type and genome filter.
    my $source = $args->{-source};
    my $genome = $args->{-genome};

    #
    # Check for the case where we use the optimized/cached version.
    #
    if ((!defined($source) || $source eq 'SEED') && !defined($genome))
    {
	return $self->_ids_to_functions_opt1($ids);
    }

    # Build the filter string, object name list, and parameter value list prefix.
    my ($objects, $filter, @parms) = $sap->ComputeFeatureFilter($source, $genome);
    # Loop through the incoming features.
    for my $id (@$ids) {
        # Get the functional assignment for this feature.
        my ($function) = $sap->GetFlat($objects, $filter, [@parms, $id],
                                       "Feature(function)");
        # Only proceed if we found one.
        if ($function) {
            Trace("Function found for $id.") if T(3);
            $retVal->{$id} = $function;
        }
    }
    # Return the result.
    return $retVal;
}

use Data::Dumper;

sub _ids_to_functions_opt1
{
    my($self, $ids) = @_;

    my $q = $self->{db}->{_dbh}->quote();
    my $out = $self->_memcache_accelerate($ids, "f", sub {
	my($self, $id_hash, $out, $upd) = @_;

	my @ids = keys %$id_hash;
	my $qs = join(", ", map { "?" } 0..$#ids);
	my $res = $self->{db}->{_dbh}->SQL(qq(SELECT id, function
					      FROM ${q}Feature$q
					      WHERE id IN ($qs)), undef, @ids);
	for my $ent (@$res)
	{
	    my($id, $fn) = @$ent;
	    $out->{$id} = $fn;
	    push(@$upd, ["f:$id", $fn, 12 * 60 * 60]) if $upd;
	}
    });
    return $out;
}

sub _memcache_accelerate
{
    my($self, $ids, $prefix, $lookup_code) = @_;

    my $mc = $self->{memcache};

    my $out = {};

    $prefix .= ":";
    my %ids = map { $_ => 1 } @$ids;
    my $update;
    if (defined($mc))
    {
	my $mcout  = $mc->get_multi(map { $prefix . $_ } keys %ids);
	# print STDERR "memcache get " . Dumper($mcout);

	for my $fid (keys %$mcout)
	{
	    my $k = $fid;
	    $fid =~ s/^$prefix//;
	    $out->{$fid} = $mcout->{$k};
	    delete $ids{$fid};
	}
	$update = [];
    }
    #
    # Look up the remaining fids
    #
    if (%ids)
    {
	my @update;
	# print STDERR "lookup " .  Dumper(\%ids);
	$lookup_code->($self, \%ids, $out, $update);
    }

    if ($update && @$update)
    {
	if ($mc->can('set_multi'))
	{
	    $mc->set_multi(@$update);
	}
	else
	{
	    $mc->set(@$_) foreach @$update;
	}
    }

    return $out;
}

sub _memcache_accelerate_list
{
    my($self, $ids, $prefix, $lookup_code) = @_;

    my $mc = $self->{memcache};

    my $out = {};

    $prefix .= ":";
    my %ids = map { $_ => 1 } @$ids;
    my $update;
    if (defined($mc))
    {
	my $mcout  = $mc->get_multi(map { $prefix . $_ } keys %ids);
	# print STDERR "memcache get " . Dumper($mcout);

	for my $fid (keys %$mcout)
	{
	    my $k = $fid;
	    $fid =~ s/^$prefix//;
	    $out->{$fid} = [split(/$;/, $mcout->{$k})];
	    delete $ids{$fid};
	}
	$update = {};
    }
    #
    # Look up the remaining fids
    #
    if (%ids)
    {
	my @update;
	# print STDERR "lookup " .  Dumper(\%ids);
	$lookup_code->($self, \%ids, $out, $update);
    }

    if ($update && %$update)
    {
	my $timeout = 12 * 60 * 60;
	if ($mc->can('set_multi'))
	{
	    my @update;
	    while (my($k, $vlist, $timeout) = each %$update)
	    {
		push(@update, [$k, join($;, @$vlist), $timeout]);
	    }

	    $mc->set_multi(@update);
	}
	else
	{
	    while (my($k, $vlist) = each %$update)
	    {
		$mc->set($k, join($;, @$vlist), $timeout);
	    }
	}
    }

    return $out;
}

=head3 occ_of_role

    my $roleHash =          $sapObject->occ_of_role({
                                -roles => [$role1, $role2, ...],
                                -functions => [$function3, $function4, ...],
                                -genomes => [$genome1, $genome2, ...],
                            });

Search for features in a specified genome with the indicated roles or
functions.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -roles (optional)

Reference to a list of the roles to search for.

=item -functions (optional)

Reference to a list of the functional assignments to search for.

=item -genomes (optional)

ID of the genomes whose genes are to be searched for the specified roles and
assignments.

=back

=item RETURN

Returns a reference to a hash that maps each specified role ID or functional
assignment to a list of the FIG IDs of genes that have that role or assignment.

    $roleHash = { $role1 => [$fid1a, $fid1b, ...],
                  $role2 => [$fid2a, $fid2b, ...],
                  $function3 => [$fid3a, $fid3b, ...],
                  $function4 => [$fid4a, $fid4b, ...],
                  ... };

=back

=cut

sub occ_of_role {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of roles.
    my $roles = ServerThing::GetIdList(-roles => $args, 1);
    # Get the list of functions.
    my $functions = ServerThing::GetIdList(-functions => $args, 1);
    # Insure we have something to look for.
    Confess("No -roles or -functions specified for occ_of_role.") if ! (@$roles + @$functions);
    # These hashes will be used to keep the results.
    my %roleHash = map { $_ => [] } @$roles;
    my %functionHash = map { $_ => {} } @$functions;
    # Get the IDs of the relevant genome.
    my $genomes = ServerThing::GetIdList(-genomes => $args, 1);
    # For backward compatability, we support the old -genome parameter.
    if (exists $args->{-genome}) {
        push @$genomes, $args->{-genome};
    }
    # Create a hash and a flag to help us filter by genome.
    my $genomeFilter = (@$genomes > 0);
    my %genomeHash = map { $_ => 1 } @$genomes;
    # We'll build the filter clause and parameters in here.
    my ($filter, @parms) = ('IsFunctionalIn(from-link) = ?');
    # If the number of genomes is small, add filtering by genome ID.
    if ($genomeFilter && @$genomes < 10) {
        # For each genome, accumulate a LIKE filter and a feature ID pattern.
        my @filters;
        for my $genome (@$genomes) {
            push @filters, 'IsFunctionalIn(to-link) LIKE ?';
            push @parms, "fig|$genome.%";
        }
        # Add the feature ID patterns to the filter string.
        $filter .= " AND (" . join(" OR ", @filters) . ")";
    }
    # We now need a list of all the roles to find in the role index. Get a hash of
    # the ones coming in directly.
    my %searchRoles = map { $_ => 1 } @$roles;
    # Add any roles inside the included functions.
    for my $function (@$functions) {
        for my $functionRole (roles_of_function($function)) {
            $searchRoles{$functionRole} = 1;
        }
    }
    # Now loop through all the roles in the search hash.
    for my $role (keys %searchRoles) {
        # Get this role's feature list.
        my $fidList = $roleHash{$role};
        # Get all of the features and their assignments for the specified role.
        my @fidData = $sap->GetAll("IsFunctionalIn Feature", $filter, [$role, @parms],
                                  [qw(Feature(id) Feature(function))]);
        # Loop through the features.
        for my $fidDatum (@fidData) {
            # Get this feature's ID and assignment.
            my ($id, $function) = @$fidDatum;
            # Check to see if it's in one of our genomes.
            my $genomeID = genome_of($id);
            if (! $genomeFilter || $genomeHash{$genomeID}) {
                # Here we want to keep this feature. If this role is in the
                # role hash, add it to the role's list.
                if (defined $fidList) {
                    push @$fidList, $id;
                }
                # If this feature's function is in the function hash, add it
                # to the function's feature hash. We're using a hash instead of
                # a list because we may see the function multiple times and need
                # to prevent duplicates.
                if (exists $functionHash{$function}) {
                    $functionHash{$function}{$id} = 1;
                }
            }
        }
    }
    # Loop through the roles, putting them in the return hash.
    for my $role (keys %roleHash) {
        $retVal->{$role} = $roleHash{$role};
    }
    # Loop through the functions, putting them in the return hash. We flatten
    # each hash of feature IDs to a list when we do this.
    for my $function (keys %functionHash) {
        $retVal->{$function} = [sort keys %{$functionHash{$function}}];
    }
    # Return the result.
    return $retVal;
}


=head2 Chemistry Methods

=head3 all_complexes

    my $complexList =       $sapObject->all_complexes();

Return a list of all the complexes in the database.

=over 4

=item RETURN

Returns a reference to a list of complex IDs.

    $complexList = [$cpx1, $cpx2, ...]

=back

=cut

sub all_complexes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of complexes.
    my @retVal = $sap->GetFlat("Complex", "", [], 'id');
    # Return the result.
    return \@retVal;
}

=head3 all_models

    my $modelHash =         $sapObject->all_models();

Return a hash of all the models in the database, mapping each one to the relevant
genome.

=over 4

=item RETURN

Returns a reference to a hash that maps each model ID to a genome ID.

    $modelHash = { $model1 => $genome1, $model2 => $genome2, ... };

=back

=cut

sub all_models {
    # Get the parameters.
    my ($self) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Read the models and genomes.
    my %retVal = map { $_->[0] => $_->[1] } $sap->GetAll('Models', "", [], "from-link to-link");
    # Return the result hash.
    return \%retVal;
}


=head3 all_reactions

    my $reactions =         $sapObject->all_reactions();

Return a list of all the reactions in the database.

=over 4

=item RETURN

Returns a reference to a list of all the reactions.

    $reactions = [$rx1, $rx2, ...];

=back

=cut

sub all_reactions {
    # Get the parameters.
    my ($self) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of reactions.
    my @retVal = $sap->GetFlat('Reaction', "", [], 'id');
    # Return it.
    return \@retVal;
}

=head3 all_roles_used_in_models

    my $rolesList =         $sapObject->all_roles_used_in_models();

Return a list of all the roles used in models.

=over 4

=item RETURN

Returns a reference to a list of role names. Each named role
triggers a complex used in at least one reaction belonging to
a model.

    $rolesList = [$role1, $role2, ...]

=back

=cut

sub all_roles_used_in_models {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get all the roles for complexes. Every complex is in at
    # least one model.
    my %roles = map { $_ => 1 } $sap->GetFlat("IsTriggeredBy",
        "", [], 'IsTriggeredBy(to-link)');
    # Return the result.
    return [keys %roles];
}

=head3 complex_data

    my $complexHash =       $sapObject->complex_data({
                                -ids => [$cpx1, $cpx2, ...],
                                -data => [$fieldA, $fieldB, ...]
                            });

Return the specified data items for each incoming reaction complex.

=over 4

=item parameter

Reference to hash with the following keys.

=over 8

=item -ids

Reference to a list of the IDs of reaction complexes of interest.

=item -data

Reference to a list of the names of the data items desired for each of the
specified complexes.

=over 12

=item name

Name of the complex (or C<undef> if the complex is nameless).

=item reactions

Reference to a list of the reactions in the complex.

=item roles

Reference to a list of 2-tuples for the roles in the complex, each
containing (0) the role name, and (1) a flag that is TRUE if the
role is optional to trigger the complex and FALSE if it is necessary.

=back

=back

=item RETURN

Returns a reference to a hash mapping each incoming complex to an n-tuple
containing the desired data fields in the order specified.

    $complexHash = { $cpx1 => [$data1A, $data1B, ...],
                     $cpx2 => [$data2A, $data2B, ...]
                     ... };

=back

=cut

sub complex_data {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Create the return hash.
    my $retVal = {};
    # Get the list of complex IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Get the list of data fields.
    my $fields = ServerThing::GetIdList(-data => $args);
    # Loop through the IDs.
    for my $id (@$ids) {
        # We'll put the data for this complex in here.
        my @data;
        # Loop through the selected data fields.
        for my $field (@$fields) {
            if ($field eq 'name') {
                # Here we need the complex name, which is taken from the Complex
                # entity.
                my ($name) = $sap->GetFlat("Complex", 'Complex(id) = ?', [$id], 'name');
                push @data, $name;
            } elsif ($field eq 'reactions') {
                # Here we need the list of reactions, which comes from the
                # IsSetOf relationship.
                my @reactions = $sap->GetFlat("IsSetOf", 'IsSetOf(from-link) = ?', [$id],
                    'to-link');
                push @data, \@reactions;
            } elsif ($field eq 'roles') {
                # The roles are taken from the IsTriggeredBy relationship.
                my @roles = $sap->GetAll("IsTriggeredBy", 'IsTriggeredBy(from-link) = ?', [$id],
                    [qw(to-link optional)]);
                push @data, \@roles;
            } else {
                # Here we have an invalid field name.
                Confess("Invalid field name $field specified in complex_data.");
            }
        }
        # Store the retrieved data in the return hash.
        $retVal->{$id} = \@data;
    }
    # Return the result.
    return $retVal;
}

=head3 coupled_reactions

    my $reactionHash =      $sapObject->coupled_reactions({
                                -ids => [$rx1, $irx2, ...]
                            });

For each of a set of reactions, get the adjacent reactions in the metabolic network.
Two reactions are considered I<adjacent> if they share at least one compound that
is neither a cofactor or a ubiquitous compound (like water or oxygen). The compounds
that relate the adjacent reactions are called the I<connecting compounds>. In most cases,
each pair of adjacent reactions will have only one connecting compound, but this is
not guaranteed to be true.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of reaction IDs.

=back

=item RETURN

Returns a reference to a hash mapping each reaction ID to a sub-hash. Each sub-hash maps
adjacent reactions to the relevant connecting compounds.

    $reactionHash = { $rx1 => { $rx1a => [$cpd1ax, $cpd1ay, ...],
                                $rx1b => [$cpd1bx, $cpd1by, ...],
                     ...};

=back

=cut

sub coupled_reactions {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of reaction IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Declare the return variable.
    my $retVal = {};
    # Loop through the IDs.
    for my $id (@$ids) {
        # This hash will be used to track the connected reactions. When we're done,
        # we'll store it in the return hash.
        my %connections;
        # Get the reactions connected to this one. For each reaction, we get the reaction ID and
        # the connecting compound.
        my @connects = $sap->GetAll('Involves Compound IsInvolvedIn',
                                    "Involves(from-link) = ? AND Compound(ubiquitous) = 0 AND Involves(cofactor) = 0 AND IsInvolvedIn(cofactor) = 0",
                                    [$id], "IsInvolvedIn(to-link) Involves(to-link)");
        # Loop through the connections.
        for my $connect (@connects) {
            # Get the adjacent reaction ID and the connection compound's ID.
            my ($reaction, $compound) = @$connect;
            # Throw away reflexive connections.
            if ($reaction ne $id) {
                # We're not reflexive, so we want to add this reaction to the connection hash.
                # Insure it isn't there already.
                if (! exists $connections{$reaction}) {
                    # List was empty, so put in this compound as a singleton.
                    $connections{$reaction} = [$compound];
                } else {
                    # The reaction is already in the hash. Pull out its list.
                    my $compoundList = $connections{$reaction};
                    # If the compound is not in there yet, add it.
                    if (! grep { $_ eq $compound } @$compoundList) {
                        push @$compoundList, $compound;
                    }
                }
            }
        }
        # Store this reaction's connection data in the return hash.
        $retVal->{$id} = \%connections;
    }
    # Return the results.
    return $retVal;
}


=head3 models_to_reactions

    my $modelHash =         $sapObject->models_to_reactions({
                                -ids => [$model1, $model2, ...]
                            });

Return the list of reactions in each specified model.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of model IDs, indicating the models of interest.

=back

=item RETURN

Returns a reference to a hash that maps each model ID to a list of the reactions in the model.

    $modelHash = { $model1 => [$rx1a, $rx1b, ...],
                   $model2 => [$rx2a, $rx2b, ...],
                   ... };

=back

=cut

sub models_to_reactions {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of model IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Declare the return hash.
    my $retVal = {};
    # Loop through the incoming model IDs.
    for my $id (@$ids) {
        # Get the list of reactions for this model.
        my @reactions = $sap->GetFlat('Requires', "Requires(from-link) = ?",
                                      [$id], 'to-link');
        # Store them in the return hash.
        $retVal->{$id} = \@reactions;
    }
    # Return the results.
    return $retVal;
}

=head3 reaction_neighbors

    my $reactionHash =      $sapObject->reactionNeighbors({
                                -ids => [$rx1, $rx2, ...],
                                -depth => 1
                            });

Return a list of the reactions in the immediate neighborhood of the specified reactions.
A separate neighborhood list will be generated for each incoming reaction; the neighborhood
will consist of reactions connected to the incoming reaction and reactions connected to those
reactions up to the specified depth. (Two reactions are I<connected> if they have a compound
in common that is not a cofactor or a ubiquitous chemical like water or ATP).

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of IDs for the reactions of interest.

=item -depth (optional)

Number of levels to which the neighborhood search should take place. If the depth is I<n>, then
the neighborhood will consist of the original reaction and every other reaction for which there is
a sequence of I<n+1> or fewer reactions starting with the original and ending with the other
reaction. Thus, if I<n> is zero, the original reaction is returned as a singleton. If I<n> is 1,
then the neighborhood is the original reaction and every reaction connected to it. The default
is C<2>.

=back

=item RETURN

Returns a reference to a hash mapping each incoming reaction to a sub-hash. The sub-hash
maps each reaction in the neighborhood to its distance from the original reaction.

    $reactionHash = { $rx1 => { $rx1a => $dist1a, $rx1b => $dist1b, ... },
                      $rx2 => { $rx2a => $dist2a, $rx2b => $dist2b, ... },
                      ... };

=back

=cut

sub reaction_neighbors {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of reaction IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Get the reaction depth. We default to 3. (Note, however, that 0 is meaningful.)
    my $depth = $args->{-depth};
    if (! defined $depth) {
        $depth = 2;
    }
    # Declare the return hash.
    my $retVal = {};
    # Loop through the incoming reaction IDs.
    for my $id (@$ids) {
        # This hash will contain the neighboring reactions found. We start with
        # the original reaction.
        my %neighbors = ($id => 0);
        # This list will contain the new reactions found in the current layer. We
        # prime the loop with the original reaction.
        my @found = ($id);
        for (my $i = 1; $i <= $depth; $i++) {
            # Get a copy of the current layer's reactions.
            my @oldFound = @found;
            @found = ();
            # Loop through the previous layer's reactions.
            for my $rx (@oldFound) {
                # Get the immediate neighbors of this reaction.
                my @neighbors = $sap->GetFlat("Involves Compound IsInvolvedIn",
                                    "Involves(from-link) = ? AND Compound(ubiquitous) = 0 AND Involves(cofactor) = 0 AND IsInvolvedIn(cofactor) = 0",
                                    [$rx], 'IsInvolvedIn(to-link)');
                # Save the new ones.
                for my $neighbor (@neighbors) {
                    if (! exists $neighbors{$neighbor}) {
                        $neighbors{$neighbor} = $i;
                        push @found, $neighbor;
                    }
                }
                Trace(scalar(@neighbors) . " neighbors found for $rx.") if T(4);
            }
            Trace(scalar(@found) . " reactions at depth $i.") if T(3);
        }
        # Save the neighbors found for this reaction.
        $retVal->{$id} = \%neighbors;
    }
    # Return the result hash.
    return $retVal;
}

=head3 reaction_path

    my $reactionList =      $sapObject->reaction_path({
                                -roles => [$role1, $role2, ...],
                                -maxLength => 10
                            });

Find the shortest reaction path that represents as many of the specified roles
as possible. Note that since the a reaction may be associated with multiple
roles, it is possible for a single role to be represented more than once in
the path.

The search is artificially limited to paths under a maximum length that can
be specified in the parameters.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -roles

Reference to a list of the roles to be covered by the reaction path.

=item -maxLength (optional)

Maximum number of reactions to allow in the reaction path. The default is two more than
the number of roles.

=back

=item RETURN

Returns a reference to a list of the best reaction paths. Each reaction path is
represented by a list of lists, the sub-lists containing the reaction IDs followed
by the roles represented by the reaction. The paths returned will be the shortest ones
found with the minimal number of missing roles.

    $reactionList = [
                     [[$rxn1a, $role1ax, $role1ay, ...], [$rxn1b, $role1bx, $role1by, ...], ...],
                     [[$rxn2a, $role2ax, $role2ay, ...], [$rxn2b, $role2bx, $role2by, ...], ...],
                     ... ];

=back

=cut

sub reaction_path {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of roles to cover.
    my $roles = ServerThing::GetIdList(-roles => $args);
    my $roleCount = scalar(@$roles);
    # Get the maximum length.
    my $maxLength = $args->{-maxLength} || $roleCount + 2;
    # To start, we need a hash that maps each reaction that relates to one of our roles
    # to a list of the related roles.
    my %rxnHash;
    for my $role (@$roles) {
        my @rxns = $sap->RoleReactions($role);
        for my $rxn (@rxns) {
            push @{$rxnHash{$rxn}}, $role;
        }
    }
    # We need to be able to create reaction paths.
    require ReactionPath;
    # Our strategy will be to build a queue of reaction paths. At any given time we will
    # track the lowest missing-roles number and the length of each path (all paths will
    # have the same length at the end of each iteration). When the missing-roles number
    # hits 0 or the length hits the maximum, we stop. This first list contains the queue
    # of paths in progress.
    my $pathList = [];
    # This hash will map each reaction we've encountered to its neighboring reactions.
    # If a reaction is not found in here, we get its information from the database.
    my %rxnNeighbors;
    # This will track the best (lowest) missing-roles number.
    my $leastMissing = $roleCount;
    # This will track the best paths found so far.
    my @bestPaths;
    # Now we prime the list with the reactions found so far.
    for my $rxn (keys %rxnHash) {
        # Get the roles covered by this reaction. There will always be at least one.
        my $foundRoles = $rxnHash{$rxn};
        # Compute the missing roles. The number of found roles is almost always one, so
        # we don't bother doing fancy hash stuff.
        my @missingRoles;
        for my $role (@$roles) {
            if (! grep { $_ eq $role } @$foundRoles) {
                push @missingRoles, $role;
            }
        }
        # Compute this path.
        my $path = ReactionPath->new($rxn, $foundRoles, \@missingRoles);
        # Compute the number of missing roles and merge them into the least-missing indicator.
        my $missing = scalar @missingRoles;
        if ($leastMissing > $missing) {
            # Here we have a new least-missing limit.
            @bestPaths = ($path);
            $leastMissing = $missing;
        } elsif ($leastMissing == $missing) {
            # Here we have a comparable path to the ones we've saved.
            push @bestPaths, $path;
        }
        # Add this path to the queue.
        push @$pathList, $path;
    }
    # Now we're ready to start the main loop. The current length of all the paths is 1.
    my $pathLength = 1;
    # This will remember the length of the shortest best-quality paths.
    my $bestPathLength = 1;
    # Loop until we hit the maximum path length or have found a path with no missing roles.
    while ($pathLength < $maxLength && $leastMissing > 0) {
        # The queue for the next iteration will be built in here.
        my @newList;
        # Update the path length. This is the length for the paths we'll be putting
        # in the new queue.
        $pathLength++;
        # Loop through the current path list.
        for my $path (@$pathList) {
            # Get the last reaction in this path.
            my $rxn0 = $path->lastReaction();
            # Look for the reactions in its neighborhood. We try to get these from the hash,
            # but if they aren't there, we query the database.
            my $neighbors = $rxnNeighbors{$rxn0};
            if (! defined $neighbors) {
                my @neighbors = $sap->GetFlat("Involves Compound IsInvolvedIn",
                                    "Involves(from-link) = ? AND Compound(ubiquitous) = 0 AND Involves(cofactor) = 0 AND IsInvolvedIn(cofactor) = 0 AND IsInvolvedIn(to-link) <> ?",
                                    [$rxn0, $rxn0], 'IsInvolvedIn(to-link)');
                $rxnNeighbors{$rxn0} = \@neighbors;
                $neighbors = \@neighbors;
            }
            # Compute the current number of missing roles for this path.
            my $currentMissing = $path->missing();
            # Loop through the neighbors, creating longer paths.
            for my $rxn1 (@$neighbors) {
                # Try to create a longer path.
                my $newPath = $path->AddReaction($rxn1, $rxnHash{$rxn1});
                if ($newPath) {
                    # Here the longer path was found. We'll set this to TRUE if we want to keep
                    # working with this path.
                    my $keep = 0;
                    # Compute the quality of this path (number of missing roles).
                    my $missing = $newPath->missing();
                    if ($leastMissing > $missing) {
                        # Here we have the new best path. Save the missing-role
                        # count and add the path to the queue.
                        $leastMissing = $missing;
                        # Re-establish the best-paths queue.
                        @bestPaths = ($newPath);
                        $bestPathLength = $pathLength;
                        # We want to keep this one.
                        $keep = 1;
                    } elsif ($leastMissing == $missing && $missing < $currentMissing &&
                             $bestPathLength == $pathLength) {
                        # Here the path is better than it used to be and is as good as
                        # and as short as the current best path. Denote it's one of the best.
                        push @bestPaths, $newPath;
                        # We want to keep it.
                        $keep = 1;
                    } elsif ($missing - $leastMissing <= $maxLength - $pathLength) {
                        # Here it's not one of the best paths, but there's a chance
                        # it can catch up if we expand to the full path length, so
                        # we keep the path, but it's not worth saving it to return to
                        # the caller.
                        $keep = 1;
                    }
                    # Put this path in the new list if we're keeping it.
                    if ($keep) {
                        push @newList, $newPath;
                    }
                }
            }
        }
        # Replace the old path list with the new one.
        $pathList = \@newList;
    }
    # We're done. Return all of the paths that have the minimal number of missing roles.
    my @retVal;
    for my $path (@bestPaths) {
        # We want to output this path. Convert the path to a list.
        my @listedPath;
        # Loop through the reactions in the path.
        for my $rxn ($path->path()) {
            # Get the roles for this reaction. If it has no roles associated with it,
            # we use an empty list.
            my $roles = $rxnHash{$rxn} || [];
            # Put the reaction and its roles into the output list for this path.
            push @listedPath, [$rxn, @$roles];
        }
        # Put the path information in the output.
        push @retVal, \@listedPath;
    }
    # Return the results found.
    return \@retVal;
}


=head3 reaction_strings

    my $reactionHash =      $sapObject->reaction_strings({
                                -ids => [$rx1, $rx2, ...],
                                -roles => 1,
                                -names => 1
                            });

Return the display string for each reaction. The display string contains the compound IDs
(as opposed to the atomic formulas) and the associated stoichiometries, with the
substrates on the left of the arrow and the products on the right.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of IDs for the reactions of interest.

=item -roles (optional)

If TRUE, then each reaction string will be associated with a list of the reaction's roles
in the result. The default is FALSE.

=item -names (optional)

If C<1>, then the compound name will be included with the ID in the output. If C<only>, the
compound name will be included instead of the ID. If C<0>, only the ID will be included. The
default is C<0>.

=back

=item RETURN

Returns a reference to a hash mapping each reaction ID to a displayable string describing the
reaction. If C<-roles> is TRUE, then instead of a string, the hash will map each reaction ID to
a list consisting of the string followed by the roles associated with the reaction.

=over 8

=item -roles FALSE

    $reactionHash = { $rx1 => $string1, $rx2 => $string2, ... }

=item -roles TRUE

    $reactionHash = { $rx1 => [$string1, $role1a, $role1b, ...],
                      $rx2 => [$string2, $role2a, $role2b, ...],
                      ...
                    }

=back

=back

=cut

sub reaction_strings {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get hte Sapling database.
    my $sap = $self->{db};
    # Get the list of reaction IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Check the options.
    my $roles = $args->{-roles} || 0;
    my $names = $args->{-names};
    # Declare the return variable.
    my $retVal = {};
    # Loop through the reactions.
    for my $id (@$ids) {
        # Get this reaction's equation.
        my @components = $sap->GetAll("Involves Compound", "Involves(from-link) = ?",
                                      [$id], "Involves(cofactor) Involves(stoichiometry) Involves(product) Compound(id) Compound(label)");
        # We'll build the substrate and product lists in here.
        my (@substrate, @product);
        # Loop through the components.
        for my $component (@components) {
            # Get the information about this component.
            my ($cofactor, $stoich, $product, $compound, $label) = @$component;
            # Compute the compound label.
            if ($names eq 'only') {
                $compound = $label;
            } elsif ($names) {
                $compound .= ": $label";
            }
            # Surround it with parentheses or brackets, depending on whether or not it's a cofactor.
            if ($cofactor) {
                $compound = "[$compound]";
            } else {
                $compound = "($compound)";
            }
            # Add the stoichiometry.
            if ($stoich != 1) {
                $compound = $stoich . $compound;
            }
            # Push the result into the appropriate list.
            if ($product) {
                push @product, $compound;
            } else {
                push @substrate, $compound;
            }
        }
        # Form the components into the result string.
        my $reactionData = join(" + ", @substrate) . " => " . join(" + ", @product);
        # Do we want roles?
        if ($roles) {
            # Get the roles and create a list out of them.
            $reactionData = [$reactionData, $sap->ReactionRoles($id)];
        }
        # Store the reaction data in the result hash.
        $retVal->{$id} = $reactionData;
    }
    # Return the result.
    return $retVal;
}

=head3 reactions_to_complexes

    my $reactionHash =      $sapObject->reactions_to_complexes({
                                -ids => [$rxn1, $rxn2, ...]
                            });

Return the complexes containing each reaction. Note that most reactions
are in more than one complex, so the complexes for each reaction are returned
as a list.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of reaction IDs for the reactions of interest.

=back

=item RETURN

Returns a reference to a hash mapping each incoming reaction to a list of the
associated complexes.

    $reactionHash = { $rxn1 => [$cpx1a, $cpx1b, ...],
                      $rxn2 => [$cpx2a, $cpx2b, ...],
                      ...
                    };

=back

=cut

sub reactions_to_complexes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the incoming ID list.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Declare the return variable.
    my $retVal = {};
    # Loop through the reactions.
    for my $id (@$ids) {
        # Get the complexes for this reaction.
        my @cpxes = $sap->GetFlat("IsElementOf",
            'IsElementOf(from-link) = ?', [$id], 'to-link');
        # Store them in the return hash.
        $retVal->{$id} = \@cpxes;
    }
    # Return the result.
    return $retVal;
}

=head3 reactions_to_roles

    my $reactionHash =      $sapObject->reactions_to_roles({
                                -ids => [$rx1, $rx2,...]
                            });

Return the roles associated with each reaction.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of reaction IDs for the reactions of interest.

=back

=item RETURN

Returns a reference to a hash mapping each incoming reaction to a list of the
associated roles.

    $reactionHash = { $rx1 => [$role1a, $role1b, ...],
                      $rx2 => [$role2a, $role2b, ...],
                      ...
                    };

=back

=cut

sub reactions_to_roles {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the incoming ID list.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Declare the return variable.
    my $retVal = {};
    # Loop through the reactions.
    for my $id (@$ids) {
        # Get the roles for this reaction.
        my @roles = $sap->ReactionRoles($id);
        # Store them in the return hash.
        $retVal->{$id} = \@roles;
    }
    # Return the result.
    return $retVal;
}

=head3 role_neighbors

    my $roleHash =          $sapObject({
                                -ids => [$role1, $role2, ...]
                            });

For each role, return a list of roles in the immediate chemical neighborhood. A role is
in the immediate chemical neighborhood of another role if the two roles are associated with
reactions that share a compound that is not ubiquitous or a cofactor.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of role names.

=back

=item RETURN

Returns a reference to a hash that maps each incoming role name to a list of the names
of the neighboring roles.

    $roleHash = { $role1 => [$role1a, $role1b, ...],
                  $role2 => [$role2a, $role2b, ...],
                  ... };

=back

=cut

sub role_neighbors {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of incoming roles.
    my $roles = ServerThing::GetIdList(-ids => $args);
    # Declare the return hash.
    my $retVal = {};
    # Loop through the roles.
    for my $role (@$roles) {
        # Get this role's neighbors. We put them in a hash to eliminate duplicates.
        my %others = map { $_ => 1 } $sap->GetFlat("Role Triggers IsSetOf Involves Compound IsInvolvedIn IsElementOf IsTriggeredBy Role2",
                        'Role(id) = ? AND Compound(ubiquitous) = 0 AND Involves(cofactor) = 0 AND IsInvolvedIn(cofactor) = 0',
                        [$role], 'Role2(id)');
        # Remove the incoming role from the result list.
        delete $others{$role};
        # Store the roles found in the result hash.
        $retVal->{$role} = [sort keys %others];
    }
    # Return the result.
    return $retVal;
}


=head3 role_reactions

    my $roleHash =          $sapObject->role_reactions({
                                -ids => [$role1, $role2, ...],
                                -formulas => 1
                            });

Return a list of all the reactions associated with each incoming role.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of role IDs for the roles of interest.

=item -formulas (optional)

If TRUE, then each reaction will be associated with its formula. The default is
FALSE, in which case for each role a simple list of reactions is returned.

=back

=item RETURN

Returns a reference to a hash, keyed by role ID. If C<-formulas> is FALSE, then
each role will map to a list of reaction IDs. If C<-formulas> is TRUE, then
each role maps to a sub-hash keyed by reaction ID. The sub-hash maps each
reaction to a chemical formula string with compound IDs in place of the
chemical labels.

=over 8

=item -formulas FALSE

    $roleHash = { $role1 => [$rxn1a, $rxn1b, ...],
                  $role2 => [$rxn2a, $rxn2b, ...},
                  ... };

=item -formulas TRUE

    $roleHash = { $role1 => { $rx1a => "$s1a1*$cpd1a1 + $s1a2*$cpd1a2 + ... => $s1ax*$cpd1ax + $s1ay*$cpd1ay + ...",
                              $rx1b => "$s1b1*$cpd1b1 + $s1b2*$cpd1b2 + ... => $s1bx*$cpd1bx + $s1by*$cpd1by + ...",
                              ... },
                  $role2 => { $rx2a => "$s2a1*$cpd2a1 + $s2a2*$cpd2a2 + ... => $s2ax*$cpd2ax + $s2ay*$cpd2ay + ...",
                              $rx2b => "$s2b1*$cpd2b1 + $s2b2*$cpd2b2 + ... => $s2bx*$cpd2bx + $s2by*$cpd2by + ...",
                              ... },
                 ... };

=back

=back

=cut

sub role_reactions {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of roles.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Find out if we want formulas.
    my $formulas = $args->{-formulas};
    # Declare the return variable.
    my $retVal = {};
    # Loop through the roles.
    for my $id (@$ids) {
        # Get the reactions for this role.
        my @reactions = $sap->RoleReactions($id);
        # Are we looking for formulas, too?
        if (! $formulas) {
            # No. Store the reactions in the return hash.
            $retVal->{$id} = \@reactions;
        } else {
            # Yes. Store the reactions and their strings.
            $retVal->{$id} = $self->reaction_strings({ -ids => \@reactions });
        }
    }
    # Return the result.
    return $retVal;
}

=head3 roles_to_complexes

    my $roleHash =          $sapObject->roles_to_complexes({
                                -ids => [$role1, $role2, ...],
                            });

Return the complexes (sets of related reactions) associated with each role in
the incoming list. Roles trigger many complexes, and a complex may be triggered
by many roles. A given role is considered either I<optional> or I<necessary>
to the complex, and an indication of this will be included in the output.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys:

=over 8

=item -ids

Reference to a list of the IDs for the roles of interest.

=back

=item RETURN

Returns a reference to a hash mapping each incoming role ID to a list of 2-tuples,
each consisting of (0) a complex ID, and (1) a flag that is TRUE if the role is
optional and FALSE if the role is necessary for the complex to trigger.

    $roleHash = { $role1 => [[$complex1a, $flag1a], [$complex1b, $flag1b], ...],
                  $role2 => [[$complex2a, $flag2a], [$complex2b, $flag2b], ...],
                  ... };

=back

=cut

sub roles_to_complexes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return hash.
    my $retVal = {};
    # Get the role IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the role IDs, getting the list of complexes for each.
    for my $id (@$ids) {
        # Get the complexes for this role.
        my @complexes = $sap->GetAll("Triggers", 'Triggers(from-link) = ?', [$id],
            [qw(to-link optional)]);
        # Store them in the return hash.
        $retVal->{$id} = \@complexes;
    }
    # Return the result hash.
    return $retVal;
}


=head2 DNA and Protein Sequence Methods

=head3 dlits_for_ids

    my $idHash =            $sapObject->dlits_for_ids({
                                -ids => [id1,id2,...],
                                -full => 1
                            });

Find the PUBMED literature references for a list of proteins. The
proteins can be specified either was FIG feature IDs or protein
sequence MD5s.

=over 4

=item parameter

The parameter should be a reference to a hash with the following
keys.

=over 8

=item -ids

Reference to a list of gene and protein IDs. For each gene,
literature references will be returned for the feature's protein.
For each protein, the literature references for the protein will
be returned. Genes should be specified using FIG feature IDs and
proteins using the MD5 of the protein sequence.

=item -full (optional)

If TRUE, then in addition to each literature article's PUBMED
ID, the article title and URL will be returned. (NOTE: these will
not always be available). The default is FALSE.

=back

=item RETURN

Returns a reference to a hash that maps each incoming ID to a list
of publications. The publications will normally be represented by
PUBMED IDs, but if C<-full> is TRUE, then each will be represented
by a 3-tuple consisting of (0) the PUBMED ID, (1) the article title,
and (2) the article URL.

=over 8

=item -full = FALSE

    $idHash = { $id1 => [$pubmed1a, $pubmed1b, ...],
                $id2 => [$pubmed2a, $pubmed2b, ...],
                ...
    };

=item -full = TRUE

    $idHash = { $id1 => [[$pubmed1a, $title1a, $url1a],
                         [$pubmed1b, $title1b, $url1b], ...],
                $id2 => [[$pubmed2a, $title2a, $url2a],
                         [$pubmed2b, $title2b, $url2b], ...],
                ...
    };

=back

=back

=cut

sub dlits_for_ids {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return hash.
    my $retVal = {};
    # Get the full-results flag.
    my $full = $args->{-full} || 0;
    # Get the list of IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through them.
    for my $id (@$ids) {
        # We'll put the object list and the filter in these variables.
        my ($objects, $filter);
        # Compute the object list and filter based on the type of ID.
        if ($id =~ /^fig/) {
            $objects = "Produces IsATopicOf";
            $filter = "Produces(from-link) = ?";
        } else {
            $objects = "IsATopicOf";
            $filter = "IsATopicOf(from-link) = ?";
        }
        # Get the PUBMED IDs for this protein.
        my @pubmeds = $sap->GetFlat($objects, $filter, [$id], "IsATopicOf(to-link)");
        # Did we find results and we're in full-results mode?
        if (@pubmeds && $full) {
            # Yes. We need to add title and link information for publications
            # where we have them. Create a filter clause to select these PUBMEDs.
            my $filter = "Publication(id) IN (" . join(", ", map { "?" } @pubmeds) .
                         ")";
            # Loop through for the publications, creating a hash of publication
            # data.
            my %pubData;
            my $q = $sap->Get("Publication", $filter, \@pubmeds);
            while (my $pub = $q->Fetch()) {
                my $pubmed = $pub->PrimaryValue("id");
                my $citation = $pub->PrimaryValue("citation");
                $pubData{$pubmed} = [$pubmed, $citation->text, $citation->link];
            }
            # Convert the pubmeds in the list to 3-tuples.
            my @newPubmeds;
            for my $pubmed (@pubmeds) {
                if ($pubData{$pubmed}) {
                    push @newPubmeds, $pubData{$pubmed};
                } else {
                    push @newPubmeds, [$pubmed, "<unknown>", "http://www.ncbi.nlm.nih.gov/pubmed/$pubmed"];
                }
            }
            # Save the modified list as our return value.
            @pubmeds = @newPubmeds;
        }
        # Store the pubmeds found in the result hash.
        $retVal->{$id} = \@pubmeds;
    }
    # Return the result.
    return $retVal;
}

=head3 equiv_ids_for_sequences

    my $labelHash =         $sapObject->equiv_ids_for_sequences({
                                -seqs => [[$label1, $comment1, $sequence1],
                                          [$label2, $comment2, $sequence2], ...]
                            });

Find all the identifiers in the database that produce the specified proteins.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -seqs

Reference to a list of protein specifications. A protein specification can be
a FASTA string, a 3-tuple consisting of (0) a label, (1) a comment,
and (2) a protein sequence, OR a 2-tuple consisting of (0) a label and (1)
a protein sequence. In other words, each specification can be a raw FASTA
string, a parsed FASTA string, or a simple [id, sequence] pair. In every case,
the protein sequence will be used to find identifiers and the label will be used
to identify the results.

=back

=item RETURN

Returns a hash mapping each incoming label to a list of identifiers from the
database that name the protein or a feature that produces the protein.

    $labelHash = { $label1 => [$id1a, $id1b, ...],
                   $label2 => [$id2a, $id2b, ...],
                   ... };

=back

=cut

sub equiv_ids_for_sequences {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the incoming sequence specifications.
    my $seqs = ServerThing::GetIdList(-seqs => $args);
    # There is the possibility the caller put in a single specifier, and we'll
    # think that it is a list of FASTA strings. To detect this, we see if
    # the first element is a string that doesn't begin with a greater-than
    # sign.
    my $first = $seqs->[0];
    if (defined $first && ! ref $first && substr($first,0,1) ne '>') {
        $seqs = [$seqs];
    }
    # Declare the return variable.
    my $retVal = {};
    # Loop through sequence specifiers.
    for my $seq (@$seqs) {
        # We need the label and the sequence for this specifier.
        my ($label, $comment, $sequence) = parse_fasta_record($seq);
        # Compute the ID for the protein in the sequence string.
        my $protID = $sap->ProteinID($sequence);
        # Find all its identifiers.
        my @ids = $sap->IdsForProtein($protID);
        # Store them in the return hash.
        $retVal->{$label} = \@ids;
    }
    # Return the results.
    return $retVal;
}

=head3 find_closest_genes

    my $nameHash =          $sapObject->find_closest_genes({
                                -genome => $genome1,
                                -seqs => { $name1 => $seq1,
                                           $name2 => #seq2,
                                           ... },
                                -protein => 1
                            });

Find the closest genes to the specified sequences in the specified genome.

Each indicated sequence will be converted to a DNA sequence and then the contigs
of the specified genome will be searched for the sequence. The genes in closest
proximity to the sequence will be returned. The sequences are named; in the return
hash, the genes found will be associated with the appropriate sequence name.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -genome

ID of the genome to search.

=item -seqs

Reference to a hash mapping names to sequences. The names will be used to associate
the genes found with the incoming sequences. DNA sequences should not contain ambiguity
characters.

=item protein (optional)

If TRUE, the sequences will be interpreted as protein sequences. If FALSE, the
sequences will be interpreted as DNA sequences.

=back

=item RETURN

Returns a reference to a hash mapping each sequence name to a list of 3-tuples, each
consisting of (0) a gene ID, (1) the location of the gene, and (2) the location of the
matching sequence.

    $nameHash = { $name1 => [[$fid1a, $loc1a, $match1a],
                             [$fid1b, $loc1b, $match1b], ...],
                  $name2 => [[$fid2a, $loc2a, $match2a],
                             [$fid2b, $loc2b, $match2b], ...],
                  ... }

=back

=cut

sub find_closest_genes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the genome ID.
    my $genome = $args->{-genome};
    Confess("No genome specified for \"find_closest_genes\".") if ! $genome;
    # Determine whether or not we have a protein sequence.
    my $protein = $args->{-protein} || 0;
    # Get the hash of sequences.
    my $seqHash = $args->{-seqs};
    Confess("No sequences specified for \"find_closest_genes\".") if ! $seqHash;
    Confess("Invalid sequence hash specified for \"find_closest_genes\".") if ref $seqHash ne 'HASH';
    # Declare the return variable.
    my $retVal = {};
    # How we process this request depends entirely on the type of search.
    if ($protein) {
        # For a protein search, we start by reading in all the pegs.
        my $pattern = "fig|$genome.peg.%";
        my %pegSequences = map { $_->[0] => $_->[1] }
            $sap->GetAll("Produces ProteinSequence",
                         'Produces(from-link) LIKE ?', [$pattern], [qw(from-link
                         ProteinSequence(sequence))]);
        # Loop through the incoming sequences.
        for my $name (keys %$seqHash) {
            # Get the sequence. We convert it to uppercase to match what's in the database.
            my $protein = uc $seqHash->{$name};
            # We'll put the hits in here.
            my @hits;
            # Loop through the pegs.
            for my $peg (keys %pegSequences) {
                # Look for the protein.
                my $pegSequence = $pegSequences{$peg};
                my $loc = 0;
                while (($loc = index($pegSequence, $protein, $loc)) >= 0) {
                    # Here we have a hit, so record it. Note that we bump
                    # the location, which both converts it to a position from
                    # an offset and insures we don't find the same subsequence
                    # again.
                    $loc++;
                    push @hits, [$peg, join(",", map { $_->String() } $sap->GetLocations($peg)),
                                       $peg . "_$loc+" . length($protein)];
                }
            }
            # If we found anything, store it in the return hash.
            if (@hits) {
                $retVal->{$name} = \@hits;
            }
        }
    } else {
        # For a DNA search, we start by reading in all the contigs.
        my @contigs = $sap->GetFlat("IsMadeUpOf", "IsMadeUpOf(from-link) = ?",
                                    [$genome], 'to-link');
        my $contigSequences = $self->contig_sequences({ -ids => \@contigs });
        # Loop through the incoming sequences.
        for my $name (keys %$seqHash) {
            # Get the sequence. We convert it to lowercase to match what's in the database
            # and fold Us to Ts.
            my $dna = lc $seqHash->{$name};
            $dna =~ tr/u/t/;
            # We'll put the hits in here.
            my @hits;
            # We need to search twice. Set up a hash to drive it.
            my %dirs = ('+' => $dna, '-' => rev_comp($dna));
            # Save its length.
            my $dnaLen = length($dna);
            # Loop through the contigs.
            for my $contig (keys %$contigSequences) {
                # Loop through the directions.
                for my $dir (keys %dirs) {
                    my $dnaSequence = $dirs{$dir};
                    Trace("Searching $contig in direction $dir for $name.") if T(3);
                    # Look for a hit in the contig.
                    my $loc = 0;
                    while (($loc = index($contigSequences->{$contig}, $dnaSequence, $loc)) >= 0) {
                        # Here we found one. Compute its location.
                        my $hitLocation = $contig . "_";
                        if ($dir eq '+') {
                            $hitLocation .= ($loc + 1) . "+";
                        } else {
                            $hitLocation .= ($loc + $dnaLen) . "-";
                        }
                        $hitLocation .= $dnaLen;
                        Trace("Hit found at $hitLocation.") if T(3);
                        # Find genes that overlap the hit region.
                        my @hitPegs = $sap->GenesInRegion($hitLocation);
                        # Loop through them, producing output.
                        for my $hitPeg (@hitPegs) {
                            # Get the peg's location.
                            my $hitPegLocation = join(",", map { $_->String() } $sap->GetLocations($hitPeg));
                            Trace("Hit peg $hitPeg is at $hitPegLocation.") if T(3);
                            # Only proceed if it goes in the correct direction.
                            if (index($hitPegLocation, $dir) >= 0) {
                                push @hits, [$hitPeg, $hitPegLocation, $hitLocation];
                            }
                        }
                        # Bump the location point to insure we don't find the same
                        # sequence again.
                        $loc++;
                    }
                }
            }
            # If any hits were found, store them in the return hash.
            if (@hits) {
                $retVal->{$name} = \@hits;
            }
        }
    }
    # Return the result hash.
    return $retVal;
}


=head3 ids_to_sequences

    my $idHash =            $sapObject->ids_to_sequences({
                                -ids => [$id1, $id2, ...],
                                -protein => 1,
                                -fasta => 1,
                                -source => 'LocusTag',
                                -genome => $genome,
                                -comments => { $id1 => $comment1,
                                               $id2 => $comment2,
                                               ... }
                            });

Compute a DNA or protein string for each incoming feature ID.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of feature IDs.

=item -source (optional)

Database source of the IDs specified-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
Use C<mixed> to allow mixed ID types (though this may cause problems when the same
ID has different meanings in different databases). Use C<prefixed> to allow IDs with
prefixing indicating the ID type (e.g. C<uni|P00934> for a UniProt ID, C<gi|135813> for
an NCBI identifier, and so forth). The default is C<SEED>.

=item -genome (optional)

ID of a specific genome. If specified, results will only be returned for genes in the
specified genome. The default is to return results for all genomes.

=item -protein (optional)

If TRUE, the output FASTA sequences will be protein sequences; otherwise, they
will be DNA sequences. The default is FALSE.

=item -fasta (optional)

If TRUE, the output sequences will be multi-line FASTA strings instead of sequences.
The default is FALSE, meaning the output sequences will be ordinary strings.

=item -comments (optional)

Allows the user to add a label or description to each FASTA formatted sequence.
The values is a reference to a hash whose keys are the ids, and the values are
the desired labels. This parameter is only used when the C<-fasta> option is
specified.

=back

=item RETURN

Returns a hash mapping the incoming IDs to sequence strings. IDs that
are not found in the database will not appear in the hash.

    $idHash = { $id1 => $sequence1, $id2 => $sequence2, ... };

=back

=cut

sub ids_to_sequences {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the ID source type.
    my $source = $args->{-source} || 'SEED';
    # Compute the ID conversion query.
    my ($idObjects, $idFilter, @idParms) = $sap->ComputeFeatureFilter($source,
                                                                      $args->{-genome});
    # Get the comment hash.
    my $comments = $args->{-comments} || {};
    # Get the strip flag. The default is stripped
    my $stripped = ($args->{-fasta} ? 0 : 1);
    # Get the sequence type.
    my $protein = $args->{-protein} || 0;
    # Extract the list of feature IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the feature IDs.
    for my $id (@$ids) {
        # Get the FIG IDs for the specified incoming ID.
        my @fids;
        if ($source eq 'SEED') {
            push @fids, $id;
        } else {
            push @fids, $sap->GetFlat($idObjects, $idFilter, [@idParms, $id],
                                      'Feature(id)');
        }
        # We'll put the sequence data we find in here.
        my @sequences;
        # Did we find any features?
        if (! @fids) {
            # No. Are we looking for proteins?
            if ($protein) {
                # Yes. Check to see if this is a protein identifier. Note that
                # in this case we're guaranteed the identifier is not a FIG ID,
                # so we can use part of the information we get from
                # ComputeFeatureFilter without worry.
                @sequences = $sap->GetFlat("Identifier Names ProteinSequence",
                                           $idFilter, [@idParms, $id],
                                           'ProteinSequence(sequence)');
            }
        }
        # Loop through the feature IDs.
        for my $fid (@fids) {
            # If the original ID and the FIG ID are different, add the
            # FIG ID to the comments.
            my $comment = $comments->{$id} || '';
            if ($id ne $fid) {
                $comment = join(" ", "[$fid]", $comment);
            }
            # Are we looking for DNA or protein data?
            if ($protein) {
                # Look for this feature's protein sequence. There should only be
                # one, so we only keep the first.
                my ($sequence) = $sap->GetFlat("Produces ProteinSequence",
                                               'Produces(from-link) = ?', [$fid],
                                               'ProteinSequence(sequence)');
                # Only proceed if we found something.
                if (defined $sequence) {
                    Trace("Creating FASTA for feature $fid.") if T(3);
                    # Create a sequence string.
                    push @sequences, create_fasta_record($id, $comment, $sequence,
                                                         $stripped);
                }
            } else {
                # We're looking for DNA, so get this feature's locations.
                my @locs = $sap->GetLocations($fid);
                # Only proceed if at least one location was found.
                if (scalar @locs) {
                    # Loop through the locations, getting DNA.
                    my $dna = join("", map { $sap->ComputeDNA($_) } @locs);
                    # Form everything into a FASTA string and store it as the result.
                    push @sequences, create_fasta_record($id, $comment, $dna, $stripped);
                }
            }
        }
        # Store the sequences in the output hash. Most of the time there's only
        # one sequence. If there's multiple, we return the first one.
        if ($stripped) {
            $retVal->{$id} = $sequences[0];
        } else {
            $retVal->{$id} = $sequences[0];
        }
    }
    # Return the result.
    return $retVal;
}

=head3 locs_to_dna

    my $locHash =           $sapObject->locs_to_dna({
                                -locations => {
                                    $label1 => $loc1,
                                    $label2 => $loc2,
                                    ... },
                                -fasta => 1
                                });

Return the DNA sequences for the specified locations.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -locations

Reference to a hash that maps IDs to locations. A location can be in the
form of a L</Location String>, a reference to a list of location strings,
a FIG feature ID, or a contig ID.

=item -fasta (optional)

If TRUE, the DNA sequences will be returned in FASTA format instead of
raw format. The default is FALSE.

=back

=item RETURN

Returns a reference to a hash that maps the incoming IDs to FASTA sequences for
the specified DNA locations. The FASTA ID will be the ID specified in the incoming
hash.

    $locHash = { $label1 => $sequence1,
                 $label2 => $sequence2,
                 ... };

=back

=cut

sub locs_to_dna {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database object.
    my $sap = $self->{db};
    # Get the ID/location map.
    my $locs = $args->{-locations};
    Confess("No locations specified.") if ! defined $locs;
    # Determine the output format.
    my $stripped = ($args->{-fasta} ? 0 : 1);
    # Loop through the incoming hash.
    for my $id (keys %$locs) {
        # Get the location string.
        my $locData = $locs->{$id};
        # We'll put our DNA in here.
        my $dna = "";
        # Determine what we have.
        if ($locData =~ /^fig\|/) {
            # Here we have a FIG ID.
            $dna = join("", map { $sap->ComputeDNA($_) } $sap->GetLocations($locData));
        } elsif (ref $locData eq 'ARRAY' || $locData =~ /^\S+_\d+[+\-_]\d+$/) {
            # Here we have a location string. We take steps to insure it is in the form
            # of a list reference.
            if (ref $locData ne 'ARRAY') {
                $locData = [$locData];
            }
            # Loop through the locations, accumulating DNA.
            for my $locString (@$locData) {
                $dna .= $sap->ComputeDNA(BasicLocation->new($locString));
            }
        } else {
            # Here we have a contig ID. Get the contig's DNA.
            my @sections = $sap->GetFlat("Contig HasSection DNASequence",
                                         'Contig(id) = ? ORDER BY DNASequence(id)', [$locData],
                                         'DNASequence(sequence)');
            $dna = join("", @sections);
        }
        # Output the DNA in FASTA form.
        $retVal->{$id} = create_fasta_record($id, undef, $dna, $stripped);
    }
    # Return the result.
    return $retVal;
}

=head3 roles_to_proteins

    my $roleHash =          $sapObject->roles_to_proteins({
                                -roles => [$role1, $role2, ...]
                            });

Return a list of the proteins associated with each of the incoming functional
roles.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -roles

Reference to a list of functional roles.

=back

=item RETURN

Returns a reference to a hash mapping each incoming role to a list of the
proteins generated by features that implement the role. The proteins will
be represented by MD5 protein IDs.

    $roleHash = { $role1 => [$prot1a, $prot1b, ...],
                  $role2 => [$prot2a, $prot2b, ...],
                  ... };

=back

=cut

sub roles_to_proteins {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of roles.
    my $roles = ServerThing::GetIdList(-roles => $args);
    # Declare the return hash.
    my $retVal = {};
    # Loop through the roles.
    for my $role (@$roles) {
        # Get the proteins for this role. We use a hash to filter out duplicates.
        my %dups;
        my @results = grep { ! $dups{$_}++ }
                      $sap->GetFlat("IsFunctionalIn Produces",
                                    'IsFunctionalIn(from-link) = ?', [$role],
                                    'Produces(to-link)');
        # Store them in the return hash.
        $retVal->{$role} = \@results;
    }
    # Return the result hash.
    return $retVal;
}

=head3 upstream

    my $featureHash =       $sapObject->upstream({
                                -ids => [$fid1, $fid2, ...],
                                -size => 200,
                                -skipGene => 1,
                                -fasta => 1,
                                -comments => { $fid1 => $comment1,
                                               $fid2 => $comment2, ...}
                            });

Return the DNA sequences for the upstream regions of the specified
features. The nucleotides inside coding regions are displayed in upper
case; others are displayed in lower case.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs of interest.

=item -size (optional)

Number of upstream nucleotides to include in the output. The default is
C<200>.

=item -skipGene (optional)

If TRUE, only the upstream region is included. Otherwise, the content
of the feature is included in the output.

=item -fasta (optional)

If TRUE, the output sequences will be multi-line FASTA strings instead of sequences.
The default is FALSE, meaning the output sequences will be ordinary strings.

=item -comments (optional)

Allows the user to add a label or description to each FASTA formatted sequence.
The values is a reference to a hash whose keys are the ids, and the values are
the desired labels. This parameter is only used when the C<-fasta> option is
specified.

=back

=item RETURN

Returns a hash mapping each incoming feature ID to the DNA sequence of
its upstream region.

    $featureHash = { $fid1 => $sequence1, $fid2 => $sequence2, ... };

=back

=cut

sub upstream {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the ID list.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Compute the options.
    my $skipGene = $args->{-skipGene} || 0;
    my $size = $args->{-size} || 200;
    my $comments = $args->{-comments} || {};
    my $stripped = ($args->{-fasta} ? 0 : 1);
    # Loop through the IDs.
    for my $fid (@$ids) {
        # Get the first location for this ID.
        my ($loc) = $sap->GetLocations($fid);
        if (defined $loc) {
            # We found a location. Get its length.
            my $locLen = $loc->Length();
            Trace("Location is " . $loc->String . " length $locLen.") if T(3);
            # Get the length of its contig.
            my $contigID = $loc->Contig;
            my $contigLen = $sap->ContigLength($contigID);
            # Extend the location over the specified upstream region.
            if ($skipGene) {
                # In skip-gene mode, we get a pure upstream location.
                $loc = $loc->Upstream($size, $contigLen);
            } else {
                # Otherwise, we simply extend upstream.
                $loc->ExtendUpstream($size, $contigLen);
            }
            Trace("Upstream location is " . $loc->String) if T(3);
            # Get the DNA for this location. It is already in lower case.
            my $dna = $sap->ComputeDNA($loc);
            Trace("DNA prefix is " . substr($dna, 0, 100) . ".") if T(3);
            # Get the direction of this location.
            my $locDir = $loc->Dir;
            # Find the other genes in this region.
            my @pegs = $sap->GenesInRegion($loc);
            Trace("Overlapping pegs are: " . join(", ", @pegs)) if T(3);
            # Get the gene locations that go in the same direction on the
            # same contig.
            my @locs;
            for my $peg (@pegs) {
                push @locs, grep { $_->Dir eq $locDir &&
                                   $_->Contig eq $contigID } $sap->GetLocations($peg);
            }
            # Loop through the locations found.
            for my $overLoc (@locs) {
                my ($start, $len) = $loc->OverlapRegion($overLoc);
                Trace("Ovelap with " . $overLoc->String . " is at $start for length $len.") if T(3);
                if ($len) {
                    $dna = substr($dna, 0, $start) . uc(substr($dna, $start, $len)) .
                           substr($dna, $start + $len);
                }
            }
            # Store the result.
            $retVal->{$fid} = create_fasta_record($fid, $comments->{$fid}, $dna, $stripped);
        }
    }
    # Return the result.
    return $retVal;
}


=head2 Expression Data Methods

=head3 all_experiments

    my $expList =           $sapObject->all_experiments();

Return a list of all the experiment names.

=over 4

=item RETURN

Returns a reference to a list of experiment names.

    $expList = [$exp1, $exp2, ...];

=back

=cut

sub all_experiments {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of experiment names.
    my @retVal = $sap->GetFlat('Experiment', '', [], 'id');
    # Return it.
    return \@retVal;
}

=head3 atomic_regulon_vectors

    my $regulonHash =       $sapObject->atomic_regulon_vectors({
                                -ids => [$ar1, $ar2, ...],
                                -raw => 0
                            });

Return a map of the expression levels for each specified atomic regulon. The
expression levels will be returned in the form of vectors with values C<-1>
(suppressed), C<1> (expressed), or C<0> (unknown) in each position. The positions
will correspond to the experiments in the order returned by L</genome_experiments>.

=over 4

=item parameter

The parameter should be a reference to a hash with the following key.

=over 8

=item -ids

Reference to a list of atomic regulon IDs.

=item -raw (optional)

If TRUE, then the vectors will be returned in the form of strings. Each string will
have the character C<+>, C<->, or space for the values 1, -1, and 0 respectively.

=back

=item RETURN

Returns a reference to a hash mapping the incoming atomic regulon IDs to the desired
vectors. The vectors will normally be references to lists of values pf 1, 0, and -1,
but they can also be represented as strings.

=over 8

=item Normal Output

    $regulonHash = { $ar1 => [$level1a, $level2a, ...],
                     $ar2 => [$level2a, $level2b, ...],
                     ... };

=item Output if -raw is TRUE

    $regulonHash = { $ar1 => $string1, $ar2 => $string2, ... };

=back

=back

=cut

sub atomic_regulon_vectors {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Decide if we're in raw mode or not.
    my $rawFlag = $args->{-raw} || 0;
    # Get the list of regulons to process.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the regulons.
    for my $id (@$ids) {
        # Get the level vectors for this regulon.
        my $qh = $sap->Get("WasGeneratedFrom", 'WasGeneratedFrom(from-link) = ? ORDER BY WasGeneratedFrom(to-link)',
                           [$id]);
        # Format them into the desired single vector.
        my $levels = ServerThing::ReadCountVector($qh, 'level-vector', $rawFlag);
        # Connect the vector to this regulon.
        $retVal->{$id} = $levels;
    }
    # Return the result.
    return $retVal;
}

=head3 atomic_regulons

    my $regulonHash =       $sapObject->atomic_regulons({
                                -id => $genome1
                            });

Return a map of the atomic regulons for the specified genome. Each atomic
regulon is a set of genes that are always regulated together. The map will
connect each regulon ID to a list of those genes. A given gene can only be
in one atomic regulon.

=over 4

=item parameter

The parameter should be a reference to a hash with the following key.

=over 8

=item -id

The ID of the genome of interest.

=back

=item RETURN

Returns a reference to a hash that maps each atomic regulon ID to a list of
the FIG IDs of its constituent genes.

    $regulonHash = { $regulon1 => [$fid1a, $fid1b, ...],
                     $regulon2 => [$fid2a, $fid2b, ...],
                     ... };

=back

=cut

sub atomic_regulons {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the genome ID.
    my $genomeID = $args->{-id};
    Confess("No genome ID specified for atomic_regulons.") if ! $genomeID;
    # Get the atomic regulon data for this genome.
    my $qh = $sap->Get("IsConfiguredBy IsFormedOf",
                       'IsConfiguredBy(from-link) = ?', [$genomeID]);
    # Loop through the result rows.
    while (my $resultRow = $qh->Fetch()) {
        # Get this regulon's ID and the ID of the feature in it.
        my $regulonID = $resultRow->PrimaryValue('IsFormedOf(from-link)');
        my $featureID = $resultRow->PrimaryValue('IsFormedOf(to-link)');
        # Put the feature into the regulon's feature list.
        push @{$retVal->{$regulonID}}, $featureID;
    }
    # Return the results.
    return $retVal;
}

=head3 coregulated_correspondence

    my $fidHash =           $sapObject->coregulated_correspondence({
                                -ids => [$fid1, $fid2, ...],
                                -pcLevel => 0.8,
                                -genomes => [$genome1, $genome2, ...]
                            });

Given a gene, return genes that may be coregulated because they correspond to
coregulated genes in genomes for which we have expression data (an
I<expression-analyzed genome>). For each incoming gene, a corresponding
gene will be found in each expression-analyzed genome. The coregulated
genes for the corresponding gene will be determined, and then these will be
mapped back to the original genome. The resulting genes can be considered
likely candidates for coregulation in the original genome.

=over 4

=item parameter

The parameter should be a reference to a hash with the following key.

=over 8

=item -ids

Reference to a list of FIG feature IDs.

=item -pcLevel (optional)

Minimum pearson coefficient level for a gene to be considered coregulated.
The default is C<0.5>.

=item -genomes (optional)

Reference to a list of genome IDs. If specified, only expression data from the
listed genomes will be used in the analysis; otherwise, all genomes with
expression data will be used.

=back

=item RETURN

Returns a reference to a hash that maps each incoming gene to a list of 4-tuples,
each 4-tuple consisting of (0) a hypothetical coregulated gene in this genome,
(1) a gene in an expression-analyzed genome corresponding to the input gene,
(2) a gene in the expression-analyzed genome coregulated with it (and that
corresponds to the hypothetical coregulated gene), and (3) the correlation
score.

    $fidHash = { $fid1 => [[$fid1a, $fid1ax, $fid1ay, $score1a],
                           [$fid1b, $fid1bx, $fid1by, $score1b],
                           ...],
                 $fid2 => [[$fid2a, $fid2ax, $fid2ay, $score2a],
                           [$fid2b, $fid2bx, $fid2by, $score2b],
                           ...],
                 ... };

=back

=cut

# Maximum number of maps to keep in memory.
use constant MAX_MAPPINGS => 500;

sub coregulated_correspondence {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the pearson coefficient level.
    my $pcLevel = $args->{-pcLevel} || 0.5;
    # Get the list of expression-data genomes.
    my $genomes;
    if (exists $args->{-genomes}) {
        $genomes = ServerThing::GetIdList(-genomes => $args);
    } else {
        $genomes = $self->expressed_genomes();
    }
    # Get the list of genes to process.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Declare the return variable.
    my $retVal = {};
    # Create a correspondence cache.
    require CorrespondenceCache;
    my $corrCache = CorrespondenceCache->new();
    # Loop through the incoming genes.
    for my $fid (@$ids) {
        # Compute the genome for this gene.
        my $idGenome = genome_of($fid);
        # The correspondence data for this gene will be put in here.
        my @coregs;
        # Loop through the expression-analyzed genomes.
        for my $xGenome (@$genomes) {
            # Get the corresponding gene in this expression-analyzed genome.
            my $xFid = $corrCache->get_correspondent($fid, $xGenome);
            # Only proceed if one was found.
            if (defined $xFid) {
                # Loop through the coregulated genes with sufficiently high pearson coefficients.
                my $qh = $sap->Get("IsCoregulatedWith",
                                   'IsCoregulatedWith(from-link) = ? AND IsCoregulatedWith(coefficient) >= ?',
                                   [$xFid, $pcLevel]);
                while (my $resultRow = $qh->Fetch()) {
                    # Get the coregulated gene in the expression-analyzed genome.
                    my $xFid2 = $resultRow->PrimaryValue('to-link');
                    # Get the pearson coefficient for the coregulation.
                    my $coefficient = $resultRow->PrimaryValue('coefficient');
                    # Look for a corresponding gene in the original genome.
                    my $fid2 = $corrCache->get_correspondent($xFid2, $idGenome);
                    # If we found one, put it in the output.
                    if (defined $fid2) {
                        push @coregs, [$fid2, $xFid, $xFid2, $coefficient];
                    }
                }
            }
        }
        # Store the hypothetically-coregulated genes found in the return hash.
        $retVal->{$fid} = \@coregs;
    }
    # Return the result.
    return $retVal;
}


=head3 coregulated_fids

    my $fidHash =           $sapObject->coregulated_fids({
                                -ids => [$fid1, $fid2, ...]
                            });

Given a gene, return the coregulated genes and their pearson coefficients.
Two genes are considered coregulated if there is some experimental evidence
that their expression levels are related: the pearson coefficient indicates
the strength of the relationship.

=over 4

=item parameter

The parameter should be a reference to a hash with the following key.

=over 8

=item -ids

Reference to a list of FIG feature IDs.

=back

=item RETURN

Returns a reference to a hash that maps each incoming FIG ID to a sub-hash.
The sub-hash in turn maps each related feature's FIG ID to its pearson
coefficient with the incoming FIG ID.

    $fidHash = { $fid1 => { $fid1a => $coeff1a, $fid1b => $coeff1b, ...},
                 $fid2 => { $fid2a => $coeff2a, $fid2b => $coeff2b, ...},
                 ... };

=back

=cut

sub coregulated_fids {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the FID IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the FIG IDs.
    for my $id (@$ids) {
        # Get the coefficient data for this FIG ID.
        my %coeffs = map { $_->[0] => $_->[1] }
                        $sap->GetAll("IsCoregulatedWith",
                                     'IsCoregulatedWith(from-link) = ?', [$id],
                                     [qw(to-link coefficient)]);
        # Put it in the result hash.
        $retVal->{$id} = \%coeffs;
    }
    # Return the results.
    return $retVal;
}

=head3 experiment_fid_levels

    my $expHash =           $sapObject->experiment_fid_levels({
                                -ids => [$exp1, $exp2, ...]
                            });

Given an experiment, return the on/off levels for all genes in that
experiment. An on/off level is either C<1> (expressed), C<-1> (inhibited),
or C<0> (unknown).

=over 4

=item parameter

The parameter should be a reference to a hash with the following key.

=over 8

=item -ids

Reference to a list of experiment IDs.

=back

=item RETURN

Returns a reference to a hash that maps each experiment ID to a sub-hash
that indicates the expression level of each gene for which the experiment
showed a result.

    $expHash = { $exp1 => { $fid1a => $level1a, $fid1b => $level1b, ... },
                 $exp2 => { $fid2a => $level2a, $fid2b => $level2b, ... },
                 ... };

=back

=cut

sub experiment_fid_levels {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the experiment IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the experiments.
    for my $id (@$ids) {
        # Create a level hash for this experiment.
        my %levels = map { $_->[0] => $_->[1] }
            $sap->GetAll("IndicatesSignalFor",
                         'IndicatesSignalFor(from-link) = ?', [$id],
                         [qw(to-link level)]);
        # Store it in the return hash.
        $retVal->{$id} = \%levels;
    }
    # Return the results.
    return $retVal;
}

=head3 experiment_regulon_levels

    my $expHash =           $sapObject->experiment_regulon_levels({
                                -ids => [$exp1, $exp2, ...]
                            });

Given an experiment, return the on/off levels for all atomic regulons
affected by that experiment. An on/off level is either C<1> (expressed), C<-1>
(inhibited), or C<0> (unknown).

=over 4

=item parameter

The parameter should be a reference to a hash with the following key.

=over 8

=item -ids

Reference to a list of experiment IDs.

=back

=item RETURN

Returns a reference to a hash that maps each experiment ID to a sub-hash
that indicates the expression level of each atomic regulon for which the
experiment showed a result.

    $expHash = { $exp1 => { $regulon1a => $level1a, $regulon1b => $level1b, ... },
                 $exp2 => { $regulon2a => $level2a, $regulon2b => $level2b, ... },
                 ... };

=back

=cut

sub experiment_regulon_levels {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the experiment IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the experiments.
    for my $id (@$ids) {
        # Create a level hash for this experiment.
        my %levels = map { $_->[0] => $_->[1] }
            $sap->GetAll("AffectsLevelOf",
                         'AffectsLevelOf(from-link) = ?', [$id],
                         [qw(to-link level)]);
        # Store it in the return hash.
        $retVal->{$id} = \%levels;
    }
    # Return the results.
    return $retVal;
}


=head3 expressed_genomes

    my $genomeList =        $sapObject->expressed_genomes((
                                -names => 1
                            });

List the IDs of genomes for which expression data exists in the database.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -names (optional)

If TRUE, then the return will be a reference to a hash mapping the genome IDs to
genome names; if FALSE, the return will be a reference to a list of genome IDs.
The default is FALSE.

=back

=item RETURN

Returns a reference to a list of genome IDs or a hash mapping genome IDs to genome names.

=over 8

=item -names FALSE

    $genomeList = [$genome1, $genome2, ...];

=item -names TRUE

    $genomeList = { $genome1 => $name1, $genome2 => $name2, ... };

=back

=back

=cut

sub expressed_genomes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get a list of all the genomes with atomic regulon data.
    my %genomes = map { $_->[0] => $_->[1] } $sap->GetAll('ProducedResultsFor Genome', '', [],
                                                'to-link Genome(scientific-name)');
    # The return value will be put in here.
    my $retVal;
    # Does the user want genome names?
    if ($args->{-names}) {
        $retVal = \%genomes;
    } else {
        $retVal = [ sort keys %genomes ];
    }
    return $retVal;
}

=head3 fid_experiments

    my $fidHash =           $sapObject->fid_experiments({
                                -ids => [$fid1, $fid2, ...],
                                -experiments => [$exp1, $exp2, ...]
                            });

Return the expression levels for the specified features in all experiments for which they
have results.

=over 4

=item parameter

The parameter should be a reference to a hash with the following key.

=over 8

=item -ids

Reference to a list of FIG feature IDs.

=item -experiments (optional)

A list of experiments. If specified, only levels from the indicated experiments will be
returned.

=back

=item RETURN

Returns a reference to a hash mapping each incoming feature ID to a list of 3-tuples,
each 3-tuple containing (0) an experiment ID, (1) the expression on/off indication (1/0/-1),
and (2) the normalized rma-value.

    $fidHash =  { $fid1 => [[$exp1a, $level1a, $rma1a],
                            [$exp1b, $level1b, $rma1b], ...],
                  $fid2 => [[$exp2a, $level2a, $rma2a],
                            [$exp2b, $level2b, $rma2b], ...],
                     ... };

=back

=cut

sub fid_experiments {
    # Get the paramters.
    my ($self, $args) = @_;
    # Get access to the Sapling database.
    my $sap = $self->{db};
    # Create the return hash.
    my $retVal = {};
    # Get the feature IDs.
    my $fids = ServerThing::GetIdList(-ids => $args);
    # Get the experiment list. This is optional, so it will return an empty list
    # if the parameter is unspecified.
    my $exps = ServerThing::GetIdList(-experiments => $args, 1);
    # Set up the experiment filter.
    my $unFiltered = (@$exps == 0);
    my %expFilter = map { $_ => 1 } @$exps;
    # Loop through the IDs.
    for my $fid (@$fids) {
        # Get the experiment list for this feature.
        my @rows = $sap->GetAll("HasIndicatedSignalFrom",
                                'HasIndicatedSignalFrom(from-link) = ?', [$fid],
                                [qw(to-link level rma-value)]);
        # Are we filtered?
        if ($unFiltered) {
            # Yes. Store it in the result hash.
            $retVal->{$fid} = \@rows;
        } else {
            # Otherwise, store a filtered list in the result hash. The experiment is in the
            # first position of the result tuple, and this is what we check against the
            # filter hash.
            $retVal->{$fid} = [ grep { $expFilter{$_->[0]} } @rows];
        }
    }
    # Return the result.
    return $retVal;
}


=head3 fid_vectors

    my $regulonHash =       $sapObject->fid_vectors({
                                -ids => [$fid1, $fid2, ...],
                                -raw => 0
                            });

Return a map of the expression levels for each specified feature (gene). The
expression levels will be returned in the form of vectors with values C<-1>
(suppressed), C<1> (expressed), or C<0> (unknown) in each position. The positions
will correspond to the experiments in the order returned by L</genome_experiments>.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs.

=item -raw (optional)

If TRUE, then the vectors will be returned in the form of strings. Each string will
have the character C<+>, C<->, or space for the values 1, -1, and 0 respectively.

=back

=item RETURN

Returns a reference to a hash mapping the incoming atomic regulon IDs to the desired
vectors. The vectors will normally be references to lists of values pf 1, 0, and -1,
but they can also be represented as strings.

=over 8

=item Normal Output

    $regulonHash = { $fid1 => [$level1a, $level2a, ...],
                     $fid2 => [$level2a, $level2b, ...],
                     ... };

=item Output if -raw is TRUE

    $regulonHash = { $fid1 => $string1, $fid2 => $string2, ... };

=back

=back

=cut

sub fid_vectors {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Decide if we're in raw mode or not.
    my $rawFlag = $args->{-raw} || 0;
    # Get the list of features to process.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the features.
    for my $id (@$ids) {
        # Get the level vectors for this regulon.
        my $qh = $sap->Get("HasLevelsFrom", 'HasLevelsFrom(from-link) = ? ORDER BY HasLevelsFrom(to-link)',
                           [$id]);
        # Format them into the desired single vector.
        my $levels = ServerThing::ReadCountVector($qh, 'level-vector', $rawFlag);
        # Connect the vector to this regulon.
        $retVal->{$id} = $levels;
    }
    # Return the result.
    return $retVal;
}

=head3 fids_expressed_in_range

    my $genomeHash =        $sapObject->fids_expressed_in_range({
                                -ids => [$genome1, $genome2, ...],
                                -minLevel => $min,
                                -maxLevel => $max
                            });

Return for each genome the genes that are expressed in a given fraction of the experiments
for that ganome.

=over 4

=item parameter

The parameter should be a reference to a hash containing the following keys.

=over 8

=item -ids

Reference to a list of IDs for the genomes of interest.

=item -minLevel (optional)

Minimum expression level. Only genes expressed at least this fraction of the time will be
output. Must be between C<0> and C<1> (inclusive) to be meaningful. The default
is C<0>, which gets everything less than or equal to the maximum level.

=item -maxLevel (optiona;)

Maximum expression level. Only genes expressed no more than this fraction of the time will be
output. Must be between C<0> and C<1> (inclusive) to be meaningful. The default is C<1>,
which gets everything greater than or equal to the minimum level.

=back

=item RETURN

Returns a hash that maps each incoming genome ID to a sub-hash. The sub-hash maps the FIG ID for
each qualifying feature to the level (as a fraction of the total experiments recorded) that it
is expressed.

    $genomeHash = { $genome1 => { $fid1a => $level1a, $fid1b => $level1b, ...},
                    $genome1 => { $fid2a => $level2a, $fid2b => $level2b, ...},
                  };

=back

=cut

sub fids_expressed_in_range {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of genome IDs.
    my $genomes = ServerThing::GetIdList(-ids => $args);
    # Get the minimum and maximum levels.
    my $minLevel = $args->{-minLevel}; $minLevel = 0.00 if ! defined $minLevel;
    my $maxLevel = $args->{-maxLevel}; $maxLevel = 1.00 if ! defined $maxLevel;
    # Declare the result hash.
    my $retVal = {};
    # Loop through the incoming genomes.
    for my $genome (@$genomes) {
        # Create the sub-hash for this genome.
        my $subHash = {};
        # Create a query to get the expression levels for this genome.
        my $qh = $sap->Get("HasLevelsFrom",
                           'HasLevelsFrom(from-link) LIKE ? ORDER BY HasLevelsFrom(from-link)',
                           ["fig|$genome%"], [qw(from-link level-vector)]);
        # There is a possibility we will get multiple results for a single gene. In that case
        # we average them together. These variables remember the previous gene, its total
        # fraction, and the number of vectors.
        my ($fid, $total, $count) = ('', 0, 0);
        # Loop through the results.
        while (my $row = $qh->Fetch()) {
            # Get the feature ID.
            my $newFid = $row->PrimaryValue('from-link');
            # Get the raw level vector. This will be a string of +, 0, and - characters.
            # We want the fraction of + characters over the total count of + and -
            # characters. If both counts are 0, then the result is ignored.
            my ($vector) = $row->Value('level-vector', 1);
            my $plusses = ($vector =~ tr/+//);
            my $minuses = ($vector =~ tr/-//);
            # Only proceed if we had at least one experiment with a definite result.
            if ($plusses > 0 || $minuses > 0) {
                my $level = $plusses / ($plusses + $minuses);
                # Is this a new feature?
                if ($newFid ne $fid) {
                    # Yes. Store the previous feature (if any).
                    if ($count > 0) {
                        my $oldLevel = $total / $count;
                        if ($oldLevel >= $minLevel && $oldLevel <= $maxLevel) {
                            $subHash->{$fid} = $oldLevel;
                        }
                    }
                    # Initialize for the new feature.
                    ($fid, $total, $count) = ($newFid, 0, 0);
                }
                # Add this data to the existing information about the feature.
                $total += $level;
                $count++;
            }
        }
        # If we have data left over, write it out to the sub-hash.
        if ($count > 0) {
            my $oldLevel = $total / $count;
            if ($oldLevel >= $minLevel && $oldLevel <= $maxLevel) {
                $subHash->{$fid} = $oldLevel;
            }
        }
        # Save this genome's results.
        $retVal->{$genome} = $subHash;
    }
    # Return the results.
    return $retVal;
}

=head3 fids_to_regulons

    my $fidHash =           $sapObject->fids_to_regulons({
                                -ids => [$fid1, $fid2, ...]
                            });

Return the atomic regulons associated with each incoming gene.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs for the genes of interest.

=back

=item RETURN

Returns a reference to a hash of hashes, keyed on FIG feature ID.
Each feature is mapped to a sub-hash that maps the feature's atomic
regulons to the number of features in each regulon.

    $fidHash = { $fid1 => { $regulon1a => $size1a, $regulon1b => $size1b, ...},
                 $fid2 => { $regulon2a => $size2a, $regulon2b => $size2b, ...},
                 ... };

=back

=cut

sub fids_to_regulons {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Declare the return hash.
    my $retVal = {};
    # Create a hash of atomic regulons. Each regulon will be mapped
    # to its size so that we don't have to ask for the same regulon
    # twice.
    my %regulons;
    # Loop through the FIG IDs.
    for my $id (@$ids) {
        # Get the regulons for this feature.
        my (@fidRegulons) = $sap->GetFlat("IsFormedInto",
            'IsFormedInto(from-link) = ?', [$id], 'to-link');
        # Loop through the regulons, insuring we know the size of
        # each one.
        for my $regulon (@fidRegulons) {
            # Get the size of this regulon if we don't already have it.
            if (! exists $regulons{$regulon}) {
                $regulons{$regulon} = $sap->GetCount("IsFormedOf",
                    "IsFormedOf(from-link) = ?", [$regulon]);
            }
        }
        # Store the regulons for this feature in the result hash.
        $retVal->{$id} = { map { $_ => $regulons{$_} } @fidRegulons };
    }
    # Return the result hash.
    return $retVal;
}


=head3 genome_experiments

    my $genomeHash =        $sapObject->genome_experiments({
                                -ids => [$genome1, $genome2, ...]
                            });

Return a list of the experiments for each indicated genome.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of genome IDs. For each genome ID, a list of relevant
experiments will be produced.

=back

=item RETURN

Returns a hash mapping each incoming genome ID to a list of experiments related
to that genome ID.

    $featureHash = { $id1 => [$exp1a, $exp1b, ...],
                     $id2 => [$exp2a, $exp2b, ...] };

=back

=cut

sub genome_experiments {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through them.
    for my $genomeID (@$ids) {
        # Get this genome's experiments.
        my @results = $sap->GetFlat("HadResultsProducedBy HasResultsIn",
                                    'HadResultsProducedBy(from-link) = ? ORDER BY HasResultsIn(from-link), HasResultsIn(sequence)',
                                    [$genomeID], 'HasResultsIn(to-link)');
        # Store them in the return hash.
        $retVal->{$genomeID} = \@results;
    }
    # Return the result.
    return $retVal;
}

=head3 genome_experiment_levels

    my $fidHash =           $sapObject->genome_experiment_levels({
                                -genome => $genome1,
                                -experiments => [$exp1, $exp2, ...]
                            });

Return the expression levels for the specified features in all experiments for which they
have results.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -genome

ID of a genome for which expression data is present.

=item -experiments (optional)

A list of experiments. If specified, only levels from the indicated experiments will be
returned.

=back

=item RETURN

Returns a reference to a hash mapping each of the genome's feature IDs to a list of 3-tuples,
each 3-tuple containing (0) an experiment ID, (1) the expression on/off indication (1/0/-1),
and (2) the normalized rma-value.

    $fidHash =  { $fid1 => [[$exp1a, $level1a, $rma1a],
                            [$exp1b, $level1b, $rma1b], ...],
                  $fid2 => [[$exp2a, $level2a, $rma2a],
                            [$exp2b, $level2b, $rma2b], ...],
                     ... };

=back

=cut

sub genome_experiment_levels {
    # Get the paramters.
    my ($self, $args) = @_;
    # Get access to the Sapling database.
    my $sap = $self->{db};
    # Create the return hash.
    my $retVal = {};
    # Get the genome ID.
    my $genome = $args->{-genome};
    Confess("No genome specified for genome_experiment_levels") if ! defined $genome;
    # Get the experiment list. This is optional, so it will return an empty list
    # if the parameter is unspecified.
    my $exps = ServerThing::GetIdList(-experiments => $args, 1);
    # Set up the experiment filter.
    my $unFiltered = (@$exps == 0);
    my %expFilter = map { $_ => 1 } @$exps;
    Trace("Retrieving experiment data for $genome.") if T(3);
    # Get the experiment data for this genome.
    my @params = ("fig|$genome.%");
    # If there is experiment filtering, add the experiment IDs to the condition.
    my $cond = "";
    if (@$exps) {
        push(@params, @$exps);
        my $qs = join(", ", map { "?" } @$exps);
        $cond = " AND HasIndicatedSignalFrom(to-link) IN ($qs)";
    }
    my @rows = $sap->GetAll("HasIndicatedSignalFrom",
                            "HasIndicatedSignalFrom(from-link) LIKE ? $cond ORDER BY HasIndicatedSignalFrom(from-link)",
                            \@params,
                            [qw(from-link to-link level rma-value)]);
    Trace(scalar(@rows) . " of experiment data found.") if T(SAP => 3);
    # Now we loop through the results, organizing them by feature ID. The current
    # feature ID will be kept in here.
    my $fid;
    # This will contain the rows for the current feature.
    my $fidRows = [];
    # Finally, we add a trailer row to insure the last feature's data gets stored
    # in the loop.
    push @rows, ["fig|TRAILER", $exps->[0]];
    for my $row (@rows) {
        # Is this row for an experiment we care about?
        if ($unFiltered || $expFilter{$row->[1]}) {
            # Yes. Split out the feature ID.
            my ($rowFid, @rowData) = @$row;
            # Is this a new feature or the same old one?
            if ($fid eq $rowFid) {
                # Same old one: queue the data.
                push @$fidRows, \@rowData;
            } else {
                # New feature. If we have data for the previous feature, write it out.
                if (@$fidRows) {
                    $retVal->{$fid} = $fidRows;
                }
                # Initialize for the new feature.
                $fidRows = [\@rowData];
                $fid = $rowFid;
            }
        }
    }
    # Return the result.
    return $retVal;
}

=head3 regulons_to_fids

    my $regHash =           $sapObject->regulons_to_fids({
                                -ids => [$regulon1, $regulon2, ...]
                            });

Return the list of genes in each specified atomic regulon.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of atomic regulon IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming atomic regulon ID
to a list of the FIG feature IDs for the genes found in the regulon.

    $regHash = { $regulon1 => [$fid1a, $fid1b, ...],
                 $regulon2 => [$fid2a, $fid2b, ...],
                 ... };

=back

=cut

sub regulons_to_fids {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return hash.
    my $retVal = {};
    # Get the list of regulon IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through them.
    for my $id (@$ids) {
        # Get the features in this regulon.
        my @fids = $sap->GetFlat("IsFormedOf", "IsFormedOf(from-link) = ?",
            [$id], "to-link");
        # Store them in the return hash.
        $retVal->{$id} = \@fids;
    }
    # Return the result hash.
    return $retVal;
}

=head2 Feature (Gene) Data Methods

NOTE: To get the functional assignment for a feature, see
L</Annotation and Assertion Data Methods>.

=head3 compared_regions

    my $result =           $sapObject->compared_regions({
                                -focus => $fid1,
                                -genomes => [$genome1, $genome2, ... ],
                                -extent => 16000
                            });

Return information about the context of a focus gene and the corresponding genes in
other genomes (known as I<pinned genes>). The information returned can be used to
create a compare-regions display.

The return information will be in the form of a reference to a list of contexts,
each context containing genes in a region surrounding the pinned gene on a particular
genome. The genome containing the focus gene will always be the first in the list.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -focus

The FIG ID of the focus gene.

=item -count (optional)

The number of pinned genes desired. If specified, the closest genes to the focus
gene will be located, at most one per genome. The default is C<4>.

=item -genomes (optional)

Reference to a list of genomes. If specified, only genes in the specified
genomes will be considered pinned.

=item -pins (optional)

Reference to a list of FIG feature IDs. The listed genes will be used as the pinned
genes. If this option is specified, it overrides C<-count> and C<-genomes>.

=item -extent (optional)

The number of base pairs to show in the context for each particular genome. The
default is C<16000>.

=back

=item RETURN

Returns a hash that maps each focus gene to the compared regions view for that gene.

Each compared regions view is a list of hashes, one hash per genome.

Each genome has the following keys:

    genome_id => this genome's id
    genome_name => this genome's name
    row_id => the row number for this genome
    features => the features for this genome.

The features lists will consist of one or more 9-tuples, one per gene in the context. Each
8-tuple will contain (0) the gene's FIG feature ID, (1) its functional assignment,
(2) its FIGfam ID, (3) the contig ID, (4) the start location, (5) the end location,
(6) the direction (C<+> or C<->), (7) the row index, and (8) the color index. All
genes with the same color have similar functions.

    $result = { focus_fid =>
	       [
		 { row_id => 0, genome_name => "g1name", genome_id => "g1id",
		   features => [[$fid1a, $function1a, $figFam1a, $contig1a, $start1a, $end1a, $dir1a, 0, $color1a],
				[$fid1b, $function1b, $figFam1b, $contig1b, $start1b, $end1b, $dir1b, 0, $color1b],
			 	... ],
		},
	        { row_id => 1, genome_name => "g2name", genome_id => "g2id",
		  features => [[$fid2a, $function2a, $figFam2a, $contig2a, $start2a, $end2a, $dir2a, 1, $color2a],
			       [$fid2b, $function2b, $figFam2b, $contig2b, $start2b, $end2b, $dir2b, 1, $color2b],
			        ... ],
		},

                ...
		]
		};

=back

=cut

sub compared_regions {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the compared-region module.
    require SapCompareRegions;
    # Get the focus gene ID.
    my $focus = $args->{-focus};
    Confess("No focus gene specified for compared_regions.") if ! $focus;
    # Now we need to determine the pinned genes. First, we check for an explicit list.
    my $pins = ServerThing::GetIdList(-pins => $args, 1);
    if (! @$pins) {
        # Here no explicit list was specified, so we have to compute the pins. We'll
        # use this argument hash to do it.
        my %pinArgs = (-focus => $focus);
        # Check for a genome list.
        my $genomes = ServerThing::GetIdList(-genomes => $args, 1);
        if (@$genomes) {
            $pinArgs{-genomes} = $genomes;
            Trace("Using genome list.") if T(SapCompareRegions => 3);
        }
        # Get the count.
        $pinArgs{-count} = $args->{-count} || 4;
        # Ask for the pins.
        $pins = SapCompareRegions::get_pin($self, \%pinArgs);
        Trace(scalar(@$pins) . " pins found.") if T(SapCompareRegions => 3);
    }
    # Now we have the pins and the focus gene. Get the extent.
    my $extent = $args->{-extent} || 16000;
    # Compute the context list.
    my $ctxList = SapCompareRegions::get_context($self, { -focus => $focus,
                                                          -pin => $pins,
                                                          -extent => $extent });
    # Add the colors.
    my $retVal = SapCompareRegions::cluster_by_function($self, {-context => $ctxList});
    $retVal = { $focus => $retVal };
    # Return the result.
    return $retVal;
}


=head3 equiv_sequence_ids

    my $idHash =            $sapObject->equiv_sequence_ids({
                                -ids => [$id1, $id2, ...],
                                -precise => 1
                            });

Return all identifiers for genes in the database that are
protein-sequence-equivalent to the specified identifiers. In this case, the
identifiers are assumed to be in their natural form (without prefixes). For
each identifier, the identified protein sequences will be found and then
for each protein sequence, all identifiers for that protein sequence or for
genes that produce that protein sequence will be returned.

Alternatively, you can ask for identifiers that are precisely equivalent, that is,
that identify the same location on the same genome.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of identifiers of interest. These can be normal feature
identifiers in prefixed form (e.g. C<cmr|NT03SD3201>, C<gi|90022544>,
C<fig|100226.1.peg.3361>) or their natural, un-prefixed form (C<NT03SD3201>,
C<90022544>). In addition, they can be protein sequence IDs formed by taking the
hexadecimal MD5 hash of the protein sequence with an optional C<md5> or
C<gnl|md5> prefix (C<500009d8cf094fa4e6a1ebb15295c60f>,
C<gnl|md5|6a00b57a9facf5056c68e5d7fe157814>).

=item -precise

If TRUE, then only identifiers that refer to the same location on the same
genome will be returned. The default is FALSE (return all sequence-equivalent
IDs). If this option is specified, identifiers that refer to proteins rather
than features will return no result.

=item -assertions

If TRUE, then instead of returning a hash of lists, this method will return
a hash of sub-hashes. Each sub-hash will be keyed by the equivalent IDs, and
will map each ID to a list of 3-tuples describing assertions about the ID,
each 3-tuple consisting of (0) an assertion of function, (1) the source of the
assertion, and (2) a flag that is TRUE for an expert assertion and FALSE
otherwise. IDs in a sub-hash which are not associated with assertions will map
to an empty list.

=back

=item RETURN

Returns a reference to a hash that maps each incoming identifier to a list
of sequence-equivalent identifiers.

=over 4

=item Normal Output

    $idHash = { $id1 => [$id1a, $id1b, ...],
                $id2 => [$id2a, $id2b, ...],
                ... };

=item Output with -assertions = 1

    $idHash = { $id1 => { $id1a => [[$assert1ax, $source1ax, $flag1ax],
                                    [$assert1ay, $source1ay, $flag1ay], ...],
                          $id1b => [[$assert1bx, $source1bx, $flag1bx],
                                    [$assert1by, $source1by, $flag1by], ...]},
                          ... },
                $id2 => { $id2a => [[$assert2ax, $source2ax, $flag2ax],
                                    [$assert2ay, $source2ay, $flag2ay], ...],
                          $id2b => [[$assert2bx, $source2bx, $flag2bx],
                                    [$assert2by, $source2by, $flag2by], ...]},
                          ... },
                ... };

=back

The output identifiers will not include protein sequence IDs: these are
allowed on input only as a convenience.

=back

=cut

sub equiv_sequence_ids {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Check for a precise-equivalence request.
    my $precise = $args->{-precise} || 0;
    # Check for an assertions request.
    my $assertions = $args->{-assertions} || 0;
    # Construct the filter clause we'll use for identifiers.
    my $filter = 'Identifier(natural-form) = ? OR Identifier(id) = ?';
    # Loop through the IDs in the list.
    for my $id (@$ids) {
        # We'll store the equivalent IDs we find in here.
        my @results;
        # Is this precise equivalence?
        if ($precise) {
            # Ask for all identifiers that connect to a feature identified by this ID.
            @results = $sap->GetFlat("Identifier Identifies Feature IsIdentifiedBy Identifier2",
                                     $filter, [$id, $id],
                                     'Identifier2(id)');
        } else {
            # We'll put the proteins of interest in here.
            my @prots;
            # Is this a protein sequence ID?
            if (my $prot = $sap->IsProteinID($id)) {
                # Use it unmodified.
                push @prots, $prot;
            } else {
                # Here we have a database ID. Ask for all the protein sequences
                # this ID identifies directly.
                push @prots, $sap->GetFlat("Identifier Names ProteinSequence", $filter,
                                           [$id, $id], 'ProteinSequence(id)');
                # Add the ones it identifies through a feature.
                push @prots, $sap->GetFlat("Identifier Identifies Feature Produces ProteinSequence",
                                           $filter, [$id, $id], 'ProteinSequence(id)');
            }
            # Loop through the protein sequences, finding equivalent IDs.
            for my $prot (@prots) {
                push @results, $sap->IdsForProtein($prot);
            }
        }
        # Loop through the results, removing duplicates.
        my %results;
        for my $result (@results) {
            $results{$result} = 1;
        }
        # Format the output according to the assertions option.
        if ($assertions) {
            for my $result (keys %results) {
                # Get the assertion data for this ID.
                my @assertRows = $sap->GetAll("Identifier HasAssertionFrom Source",
                                              'Identifier(id) = ? ',
                                              [$result],
                                              [qw(HasAssertionFrom(function)
                                                  Source(id)
                                                  HasAssertionFrom(expert))]);
                # Store them in the hash.
                $results{$result} = \@assertRows;
            }
            # Attach the hash to this ID.
            $retVal->{$id} = \%results;
        } else {
            # Normal mode. Store the IDs found as a list.
            $retVal->{$id} = [ sort keys %results ];
        }
    }
    # Return the result.
    return $retVal;
}

=head3 fid_correspondences

    my $featureHash =       $sapObject->fid_correspondences({
                                -ids => [$fid1, $fid2, ...],
                                -genomes => [$genome1, $genome2, ...]
                            });

Return the corresponding genes for the specified features in the specified genomes.
The correspondences are determined in the same way as used by L</gene_correspondence_map>,
but this method returns substantially less data.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs.

=item -genomes

Reference to a list of genome IDs. For each incoming feature ID, the corresponding
features in the specified genomes will be returned.

=back

=item RETURN

Returns a reference to a hash that maps each incoming feature ID to a list of
corresponding feature IDs in the specified genomes. If no sufficiently corresponding
feature is found in any of the genomes, the feature ID will map to an empty list.

    $featureHash = { $fid1 => [$fid1a, $fid1b, ...],
                     $fid2 => [$fid2a, $fid2b, ...],
                     ... };

=back

=cut

sub fid_correspondences {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get feature and genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    my $genomes = ServerThing::GetIdList(-genomes => $args);
    # Sort the incoming IDs by genome.
    my %idGroups;
    for my $id (@$ids) {
        my $genome = genome_of($id);
        push @{$idGroups{$genome}}, $id;
    }
    # Loop through the incoming IDs, one genome at a time.
    for my $sourceGenome (keys %idGroups) {
        # Get the features for this source genome.
        my $sourceFids = $idGroups{$sourceGenome};
        # Loop through the target genomes.
        for my $targetGenome (@$genomes) {
            # We only need to do this if the source and target genomes are different.
            if ($sourceGenome ne $targetGenome) {
                # Get the correspondences for this genome.
                my $geneHash = $self->gene_correspondence_map({ -genome1 => $sourceGenome,
                                                                -genome2 => $targetGenome });
                # Put all the useful results into the return hash.
                for my $id (@$sourceFids) {
                    if (exists $geneHash->{$id}) {
                        push @{$retVal->{$id}}, $geneHash->{$id};
                    }
                }
            }
        }
    }
    # Return the result hash.
    return $retVal;
}


=head3 fid_locations

    my $featureHash =       $sapObject->fid_locations({
                                -ids => [$fid1, $fid2, ...],
                                -boundaries => 1
                            });

Return the DNA locations for the specified features.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs.

=item -boundaries (optional)

If TRUE, then for any multi-location feature, a single location encompassing all the
location segments will be returned instead of a list of all the segments. If the
segments cross between contigs, then the behavior in this mode is undefined
(something will come back, but it may not be what you're expecting). The default
is FALSE, in which case the locations for each feature will be presented in a list.

=back

=item RETURN

Returns a reference to a hash mapping each feature ID to a list of location strings
representing the feature locations in sequence order.

    $featureHash = { $fid1 => [$loc1a, $loc1b, ...],
                     $fid2 => [$loc2a, $loc2b, ...],
                     ... };

=back

=cut

sub fid_locations {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database object.
    my $sap = $self->{db};
    # Determine the operating mode (boundaries or normal).
    my $boundaryMode = $args->{-boundaries} || 0;
    # Get the list of identifiers.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through them.
    for my $fid (@$ids) {
        # Get the list of locations for this feature.
        my @locs = $sap->GetLocations($fid);
        # Convert the locations to location strings.
        my @locStrings = map { $_->String() } @locs;
        # Only proceed if we found something.
        if (scalar @locs) {
            # Process according to the output mode.
            if ($boundaryMode) {
                # Here we're looking for boundaries.
                my ($contig, $min, $max) = boundaries_of(\@locStrings);
                # Get the first location's direction.
                my $dir = $locs[0]->Dir;
                # Compute a location from the boundaries.
                my $locLen = $max - $min + 1;
                my $retLoc = $contig . "_";
                if ($dir eq '+') {
                    $retLoc .= "$min+$locLen";
                } else {
                    $retLoc .= "$max-$locLen";
                }
                # Store it in the return hash.
                $retVal->{$fid} = $retLoc;
            } else {
                # Here we want the list of locations, a much simpler operation.
                $retVal->{$fid} = \@locStrings;
            }
        }
    }
    # Return the result.
    return $retVal;
}

=head3 fid_map_for_genome

    my $idHash =            $sapObject->get_map_for_genome({
                                -idHash => { $myID1 => [$id1a, $id1b, ...],
                                             $myID2 => [$id2a, $id2b, ...],
                                             ... },
                                -genome => $genome1
                            });

Find FIG IDs corresponding to caller-provided genes in a specific genome.

In some situations you may have multiple external identifiers for
various genes in a genome without knowing which ones are present in the Sapling
database and which are not. The external identifiers present in the Sapling
database are culled from numerous sources, but different genomes will tend to
have coverage from different identifier types: some genomes are represented
heavily by CMR identifiers and have no Locus Tags, others have lots of Locus
Tags but no CMR identifiers, and so forth. This method allows you to throw everything
you have at the database in hopes of finding a match.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -idHash

Reference to a hash that maps caller-specified identifiers to lists of external
identifiers in prefixed form (e.g. C<LocusTag:SO1103>, C<uni|QX8I1>, C<gi|4808340>).
Each external identifier should be an alternate name for the same gene.

=item -genome (optional)

ID of a target genome. If specified, only genes in the specified target genome
will be returned.

=back

=item RETURN

Returns a hash mapping the original caller-specified identifiers to FIG IDs in the
target genome. If the identifier list is ambiguous, the first matching FIG ID will
be used. If no matching FIG ID is found, an undefined value will be used.

    $idHash = { $myID1 => $fid1, $myID2 => $fid2, ... };

=back

=cut

sub fid_map_for_genome {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the genome ID.
    my $genome = $args->{-genome};
    # Get the ID hash.
    my $idHash = $args->{-idHash};
    if (! defined $idHash) {
        Confess("No id hash specified to fid_map_for_genome.");
    } elsif (ref $idHash ne 'HASH') {
        Confess("Invalid id hash specified to fid_map_for_genome.");
    } else {
        # Here we have a valid ID hash. First, compute the object list and filter
        # clause for the ID query.
        my ($idObjects, $idFilter, @idParms) = $sap->ComputeFeatureFilter('prefixed',
                                                                          $genome);
        # Loop through the keys of the hash. These are the user's preferred IDs.
        for my $id (keys %$idHash) {
            # Get the identifiers for this ID. Note we are prepared for a singleton
            # instead of a list.
            my $idList = $idHash->{$id};
            $idList = [$idList] if ref $idList ne 'ARRAY';
            # Try to find a FIG ID for each identifier.
            my $foundFid;
            for my $external (@$idList) { last if $foundFid;
                # Check this external for a matching FIG ID.
                ($foundFid) = $sap->GetFlat($idObjects, $idFilter, [@idParms, $external],
                                            'Feature(id)');
            }
            # Store the FIG ID found (if any).
            $retVal->{$id} = $foundFid;
        }
    }
    # Return the result.
    return $retVal;
}


=head3 fid_possibly_truncated

    my $featureHash =       $sapObject->fid_possibly_truncated({
                                -ids => [$fid1, $fid2, ...],
                                -limit => 300
                            });

For each specified gene, return C<stop> if its end is possibly truncated,
C<start> if its beginning is possibly truncated, and an empty string
otherwise. Truncation occurs if the gene is located near either edge of a
contig.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG gene IDs.

=item -limit (optional)

The distance from the end of a contig considered to be at risk for truncation.
the default is 300.

=back

=item RETURN

Returns a hash mapping each incoming gene ID to the appropriate value (C<start>
if it has a possibly-truncated start, C<stop> if it has a possibly-truncated stop,
or the empty string otherwise). Note that the empty string is expected to be the
most common result.

    $featureHash = { $fid1 => $note1, $fid2 => $note2, ... };

=back

=cut

sub fid_possibly_truncated {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the truncation limit.
    my $limit = $args->{-limit} || 300;
    # Get the list of feature IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through them.
    for my $fid (@$ids) {
        # Get this feature's location list.
        my @locs = $sap->GetLocations($fid);
        # Compute the boundaries.
        my ($contig, $min, $max) = boundaries_of([map { $_->String } @locs]);
        # Only proceed if the location was valid.
        if (defined $contig) {
            # It's okay, so get the direction.
            my $dir = $locs[0]->Dir;
            # Get the contig length.
            my ($contigLen) = $sap->GetEntityValues(Contig => $contig, ['length']);
            Confess("Database error: contig $contig not found.") if ! defined $contigLen;
            # Determine whether we're near one of the ends.
            my $nearLeft = $min < $limit;
            my $nearRight = $max > $contigLen - $limit;
            # Compute the truncation indicator. Note that STOP truncation has priority
            # over START truncation.
            my $truncated = '';
            if ($nearLeft && $dir eq '-' || $nearRight && $dir eq '+') {
                $truncated = 'stop';
            } elsif ($nearLeft && $dir eq '+' || $nearRight && $dir eq '-') {
                $truncated = 'start';
            }
            # Store the indicator.
            $retVal->{$fid} = $truncated;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 fids_to_ids

    my $featureHash =       $sapObject->fids_to_ids({
                                -ids => [$fid1, $fid2, ...],
                                -types => [$typeA, $typeB, ...],
                                -protein => 1
                            });

Find all aliases and/or synonyms for the specified FIG IDs. For each FIG
ID, a hash will be returned that maps each ID type to a list of the IDs
of that type.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs of interest,

=item -types (optional)

Reference to a list of permissible ID types. Only ID types in this list will
be present in the output. If omitted, all ID types are permissible.

=item -protein (optional)

If TRUE, then IDs for features with equivalent protein sequences will be
returned; otherwise, only IDs for precisely equivalent genes will be returned.
The default is FALSE

=item -natural (optional)

If TRUE, then the IDs will be returned in their natural form; otherwise, the
IDs are returned in prefixed form. The default is FALSE.

=back

=item RETURN

Returns a reference to a hash that maps each feature ID to a sub-hash. Each
sub-hash maps an ID type to a list of equivalent IDs of that type.

    $featureHash = { $fid1 => { $typeA => [$id1A1, $id1A2, ...],
                                $typeB => [$id1B1, $id1B2, ...],
                                ... },
                     $fid2 =>  { $typeA => [$id2A1, $id2A2, ...],
                                 $typeB => [$id2B1, $id2B2, ...],
                                ... },
                     ... };

=back

=cut

sub fids_to_ids {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the protein mode flag.
    my $protein = ($args->{-protein} ? 1 : 0);
    # Find out if we're natural or prefixed.
    my $idField = ($args->{-natural} ? "Identifier(natural-form)" : "Identifier(id)");
    # Create the list of permissible ID types. We will set $restrictedTypes
    # to TRUE if type restrictions matter.
    my ($restrictedTypes, %permissibleTypes);
    if (defined $args->{-types}) {
        $restrictedTypes = 1;
        my $types = $args->{-types};
        if (! ref $types) {
            # A scalar was specified, so we treat it as the only allowed type.
            $permissibleTypes{$types} = 1;
        } elsif (ref $types eq 'ARRAY'){
            # Here we have an array reference, which is what we are expecting.
            %permissibleTypes = map { $_ => 1 } @$types;
        } else {
            # This is an error.
            Confess("Invalid \"types\" parameter specified: it must be a scalar or a list, but found " .
                    ref $types . " instead.");
        }
    }
    # Finally, we get the the ID list.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs.
    for my $fid (@$ids) {
        Trace("Processing $fid.") if T(3);
        # Get the aliases for this ID. We'll put them in this hash, with each
        # alias mapped to its type.
        my %aliasHash;
        # First, we get the list of exact matches.
        my @idPairs = $sap->GetAll('IsIdentifiedBy Identifier',
                                   'IsIdentifiedBy(from-link) = ?',
                                   [$fid], "$idField Identifier(source)");
        # If we are doing protein equivalence, get the protein aliases as well.
        if ($protein && (my $protID = $sap->IdentifiedProtein($fid))) {
            # We have a protein ID, so we want to get all the protein identifiers
            # AND all the feature identifiers for features with this protein.
            push @idPairs, $sap->GetAll('IsNamedBy Identifier',
                                        'IsNamedBy(from-link) = ?', [$protID],
                                        "$idField Identifier(source)");
            push @idPairs, $sap->GetAll('IsProteinFor IsIdentifiedBy Identifier',
                                        'IsProteinFor(from-link) = ?', [$protID],
                                        "$idField Identifier(source)");
        }
        # Only proceed if we found something.
        my $count = scalar @idPairs;
        Trace("$count identifiers found for $fid with protein flag $protein.") if T(3);
        if ($count) {
            # Now we have a list of (identifier, source) pairs for this feature.
            # We want to convert this into a hash of lists, mapping source types to
            # identifiers from that source. This hash is used to prevent duplicates.
            my %idsFound;
            # This hash will contain the ID lists.
            my %idLists;
            # Loop through the pairs.
            for my $idPair (@idPairs) {
                # Get the ID and its source type.
                my ($id, $source) = @$idPair;
                # Only proceed if this ID is new.
                if (! exists $idsFound{$id}) {
                    $idsFound{$id} = 1;
                    # Do we want to keep IDs of this type?
                    if (! $restrictedTypes || $permissibleTypes{$source}) {
                        # Yes. Put it in the type's list.
                        push @{$idLists{$source}}, $id;
                    }
                }
            }
            # Put our hash of IDs in the return value.
            $retVal->{$fid} = \%idLists;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 fids_to_proteins

    my $fidHash =           $sapObject->fids_to_proteins({
                                -ids => [$fid1, $fid2, ...],
                                -sequence => 1
                            });

Return the ID or amino acid sequence associated with each specified gene's protein. If the gene
does not produce a protein, it will not be included in the output.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs, representing the features of interest.

=item -sequence (optional)

If TRUE, then the output will include protein sequences; otherwise, the output will include
MD5 protein IDs. The default is FALSE.

=back

=item RETURN

Returns a reference to a hash keyed by feature ID. If C<-sequence> is FALSE, then the hash
maps each feature ID to the MD5 ID of the relevant gene's protein sequence. If C<-sequence>
is TRUE, then the hash maps each feature ID to the relevant protein sequence itself.

=over 8

=item -sequence TRUE

    $fidHash = { $fid1 => $sequence1, $fid2 => $sequence2, ... };

=item -sequence FALSE

    $fidHash = { $fid1 => $md5id1, $fid2 => $md5id2, ... };

=back

=back

=cut

sub fids_to_proteins {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of feature IDs.
    my $fids = ServerThing::GetIdList(-ids => $args);

    # Find out if we're going for sequences or IDs. This determines our output fields and
    # the query path.

    if ($self->{memcache})
    {
	if (!$args->{-sequence})
	{
	    return $self->_fids_to_proteins_opt1($fids);
	}
    }

    my ($path, $field);
    if ($args->{-sequence}) {
        # Here we want the protein sequences themselves.
        $path = "Produces ProteinSequence";
        $field = 'ProteinSequence(sequence)';
    } else {
        # Here we only care about the protein ID.
        $path = "Produces";
        $field = 'to-link';
    }
    # Declare the return hash.
    my $retVal = {};
    # Loop through the incoming IDs.
    for my $fid (@$fids) {
        # Get the protein data for this feature. There can be at most one result.
        my ($result) = $sap->GetFlat($path, "Produces(from-link) = ?", [$fid], $field);
        # Store it in the return hash.
        if ($result) {
            $retVal->{$fid} = $result;
        }
    }
    # Return the result.
    return $retVal;
}

sub _fids_to_proteins_opt1
{
    my($self, $ids) = @_;

    my $out = $self->_memcache_accelerate($ids, "f2md5", sub {
	my($self, $id_hash, $out, $upd) = @_;

	my @ids = keys %$id_hash;
	my $qs = join(", ", map { "?" } 0..$#ids);
	my $res = $self->{db}->{_dbh}->SQL(qq(SELECT to_link, from_link
					      FROM IsProteinFor
					      WHERE to_link IN ($qs)), undef, @ids);
	for my $ent (@$res)
	{
	    my($id, $md5) = @$ent;
	    $out->{$id} = $md5;
	    push(@$upd, ["f2md5:$id", $md5, 12 * 60 * 60]) if $upd;
	}

    });
    return $out;
}

=head3 fids_with_evidence_codes

    my $featureHash =       $sapObject->fids_with_evidence_codes({
                                -codes => [$code1, $code2, ...],
                                -genomes => [$genome1, $genome2, ...]
                            });

Return the ID, assignment, and evidence for all features having an
evidence code of one of the specified types. The output can be restricted
to one or more specified genomes.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -codes

Reference to a list of evidence code types. This is only the prefix, not a
full-blown code. So, for example, C<ilit> would be used for indirect literature
references, C<dlit> for direct literature references, and so forth.

=item -genomes (optional)

Reference to a list of genome IDs. If no genome IDs are specified, all features
in all genomes will be processed.

=back

=item RETURN

Returns a hash mapping each feature to a list containing the function followed by
all of the feature's evidence codes.

    $featureHash = { $fid1 => [$function1, $code1A, $code1B, ...],
                     $fid2 => [$function2, $code2A, $code2B, ...],
                     ... };

=back

=cut

sub fids_with_evidence_codes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Determine the genome list.
    my $genomes = $args->{-genomes};
    if (! defined $genomes) {
        # No genomes we specified, so we do them all.
        $genomes = [ $sap->GetFlat("Genome", "", [], "id") ];
    } elsif (! ref $genomes) {
        # A scalar genome ID is converted to a list for convenience.
        $genomes = [$genomes];
    }
    Trace(scalar(@$genomes) . " genomes selected.") if T(3);
    # Get the evidence code list.
    my $codes = ServerThing::GetIdList(-codes => $args);
    Trace(scalar(@$codes) . " evidence code types selected.") if T(3);
    # Declare the return variable.
    my $retVal = {};
    # Loop through the genomes.
    for my $genomeID (@$genomes) {
        # Loop through the evidence code.
        for my $code (@$codes) {
            Trace("Processing $code for $genomeID.") if T(3);
            # Query the database for this genome and code.
            my $qh = $sap->Get("Genome IsOwnerOf Feature",
                               'Genome(id) = ? AND Feature(evidence-code) LIKE ?',
                               [$genomeID, "$code%"]);
            # Loop through the results.
            while (my $resultRow = $qh->Fetch()) {
                # Get the data for this feature.
                my $featureId = $resultRow->PrimaryValue('Feature(id)');
                my $featureFunction = $resultRow->PrimaryValue('Feature(function)');
                my @featureEvidenceCodes = $resultRow->Value('Feature(evidence-code)');
                # Put the data in the return hash.
                $retVal->{$featureId} = [$featureFunction, @featureEvidenceCodes];
            }
        }
    }
    # Return the result.
    return $retVal;
}

=head3 genes_in_region

    my $locHash =           $sapObject->genes_in_region({
                                -locations => [$loc1, $loc2, ...],
                                -includeLocation => 1
                            });

Return a list of the IDs for the features that overlap the specified
regions on a contig.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -locations

Reference to a list of location strings (e.g. C<360108.3:NZ_AANK01000002_264528_264007>
or C<100226.1:NC_003888_3766170+612>). A location string consists of a contig ID
(which includes the genome ID), an underscore, a begin offset, and either an underscore
followed by an end offset or a direction (C<+> or C<->) followed by a length.

=item -includeLocation

If TRUE, then instead of mapping each location to a list of IDs, the hash will map
each location to a hash reference that maps the IDs to their locations.

=back

=item RETURN

Returns a reference to a hash mapping each incoming location string
to a list of the IDs for the features that overlap that location.

    $locHash = { $loc1 => [$fid1A, $fid1B, ...],
                 $loc2 => [$fid2A, $fid2B, ...],
                 ... };

=back

=cut

sub genes_in_region {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Check for the includeLocation option.
    my $includeLocation = $args->{-includeLocation} || 0;
    # Get the list of location strings.
    my $locs = ServerThing::GetIdList(-locations => $args);
    # Loop through the locations.
    for my $loc (@$locs) {
        # Get the genes in the region.
        my @fids = $sap->GenesInRegion($loc);
        # If this is include-location mode, add the location strings.
        if ($includeLocation) {
            # Loop through the features found. We'll put the data we find in here.
            my %fidData;
            for my $fid (@fids) {
                $fidData{$fid} = [ map { $_->String() } $sap->GetLocations($fid) ];
            }
            # Store the list of lists in the output hash.
            $retVal->{$loc} = \%fidData;
        } else {
            # In normal mode, the feature list goes in the output hash unaltered.
            $retVal->{$loc} = \@fids;
        }
    }
    # Return the result.
    return $retVal;
}


=head3 ids_to_data

    my $featureHash =       $sapObject->ids_to_data({
                                -ids => [$id1, $id2, ...],
                                -data => [$fieldA, $fieldB, ...],
                                -source => 'UniProt'
                            });

Return the specified data items for the specified features.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of gene identifiers. Normally, these would be FIG
feature IDs, but other identifier types can be specified if you use the
C<-source> option.

=item -data

Reference to a list of data field names. The possible data field names are
given below.

=over 12

=item evidence

Comma-delimited list of evidence codes indicating the reason for the gene's
current assignment.

=item fig-id

The FIG ID of the gene.

=item function

Current functional assignment.

=item genome-name

Name of the genome containing the gene.

=item length

Number of base pairs in the gene.

=item location

Comma-delimited list of location strings indicated the location of the gene
in the genome. A location string consists of a contig ID, an underscore, the
starting offset, the strand (C<+> or C<->), and the number of base pairs.

=item publications

Comma-delimited list of PUBMED IDs for publications related to the gene.

=back

=item -source (optional)

Database source of the IDs specified-- e.g. C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
Use C<mixed> to allow mixed ID types (though this may cause problems when the same
ID has different meanings in different databases). Use C<prefixed> to allow IDs with
prefixing indicating the ID type (e.g. C<uni|P00934> for a UniProt ID, C<gi|135813> for
an NCBI identifier, and so forth). The default is C<SEED>.

=item -genome (optional)

ID of a specific genome. If specified, results will only be returned for genes in the
specified genome. The default is to return results for genes in all genomes.

=back

=item RETURN

Returns a hash mapping each incoming ID to a list of tuples, There will be one
tuple for each feature identified by the incoming ID (because some IDs are
ambiguous there may be more than one), and the tuple will contain the
specified data fields for the computed gene in the specified order.

    $featureHash = { $id1 => [$tuple1A, $tuple1B, ...],
                     $id2 => [$tuple2A, $tuple2B, ...],
                     ... };

=back

=cut

sub ids_to_data {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # This hash is used to cache genome names for performance.
    my %genomes;
    # Get the list of fields and perform a basic validation so we know we have a
    # list of data items.
    my $fields = $args->{-data};
    Confess("No data fields specified in \"ids_to_data\".") if ! defined $fields;
    Confess("Invalid data field list in \"ids_to_data\".") if ref $fields ne 'ARRAY';
    # Get the list of feature IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Create the feature filter for IDs of this type.
    my ($objects, $filter, @parms) = $sap->ComputeFeatureFilter($args->{-source},
                                                                $args->{-genome});
    # Loop through the identifiers.
    for my $id (@$ids) {
        # The output tuples for this identifier will be put in here.
        my @tuples;
        # Get the features for this identifier.
        my @dbObjects = $sap->GetList($objects, $filter, [@parms, $id]);
        # Loop through the features found.
        for my $feature (@dbObjects) {
            # The current tuple will be built in here.
            my @tuple;
            # Get the current feature ID.
            my $fid = $feature->PrimaryValue('Feature(id)');
            # Loop through the fields we need.
            for my $field (@$fields) {
                if ($field eq 'evidence') {
                    # We do a join here because there may be multiple evidence codes.
                    push @tuple, join(", ", $feature->Value('Feature(evidence-code)'));
                } elsif ($field eq 'fig-id') {
                    # The FIG ID was extracted above.
                    push @tuple, $fid;
                } elsif ($field eq 'function') {
                    # The assignment is a field in the Feature record.
                    push @tuple, $feature->Value('Feature(function)');
                } elsif ($field eq 'genome-name') {
                    # For genome names, we need to parse the feature ID.
                    my $genomeID = genome_of($fid);
                    # If we already have this genome's name, we reuse it;
                    # otherwise, we query the database.
                    if (! $genomes{$genomeID}) {
                        ($genomes{$genomeID}) = $sap->GetEntityValues(Genome => $genomeID,
                                                                      ['scientific-name']);
                    }
                    push @tuple, $genomes{$genomeID};
                } elsif ($field eq 'length') {
                    # This is the sequence-length field from the feature record.
                    push @tuple, $feature->Value('Feature(sequence-length)');
                } elsif ($field eq 'location') {
                    # Sapling has a custom method for getting locations.
                    my @locs = $sap->GetLocations($fid);
                    # We translates the location objects to location strings and
                    # join them with commas before returning them. Note that in
                    # most cases, however, there will only be one.
                    push @tuple, join(", ", map { $_->String } @locs);
                } elsif ($field eq 'publications') {
                    # The publication data is kept in the evidence codes. For each
                    # publication, there will be a "dlit" evidence code relating to
                    # it. Immediately after the "dlit" will be a PUBMED number
                    # enclosed in parentheses. Note that in array context, the match
                    # operator obligingly returns an empty list if a match fails and
                    # a list of the parenthesized matched text if the match works.
                    push @tuple, join(", ", map { $_ =~ /dlit\((\d+)/ } $feature->Value('Feature(evidence-code)'));
                } else {
                    Confess("Invalid data field name \"$field\" in \"ids_to_data\".");
                }
            }
            # Add this feature's data to the output list for this ID.
            push @tuples, \@tuple;
        }
        # Store this ID's results.
        $retVal->{$id} = \@tuples;
    }
    # Return the result.
    return $retVal;
}

=head3 ids_to_fids

    my $idHash =            $sapObject->ids_to_fids({
                                -ids => [$id1, $id2, ...],
                                -protein => 1,
                                -genomeName => $genusSpeciesString,
                                -source => 'UniProt'
                            });

Return a list of the FIG IDs corresponding to each of the specified
identifiers. The correspondence can either be gene-based (same feature)
or sequence-based (same protein).

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of identifiers.

=item -source

Database source of the IDs specified-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
Use C<mixed> to allow mixed ID types (though this may cause problems when the same
ID has different meanings in different databases). Use C<prefixed> to allow IDs with
prefixing indicating the ID type (e.g. C<uni|P00934> for a UniProt ID, C<gi|135813> for
an NCBI identifier, and so forth).

=item -protein (optional)

If TRUE, then all FIG IDs for equivalent proteins will be returned. The default is
FALSE, meaning that only FIG IDs for the same gene will be returned.

=item -genomeName (optional)

The full or partial name of a genome or a comma-delimited list of genome IDs.
This parameter is useful for narrowing the results when a protein match is
specified. If it is omitted, no genome filtering is performed.

=back

=item RETURN

Returns a reference to a hash mapping each incoming identifier to a list
of equivalent FIG IDs.

    $idHash = { $id1 => [$fid1A, $fid1B, ...],
                $id2 => [$fid2A, $fid2B, ...],
                ... };

=back

=cut

sub ids_to_fids {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Get the source.
    my $source = $args->{-source};
    Confess("No -source specified on \"ids_to_fids\".") if ! defined $source;
    # Compute the feature filter.
    my ($objects, $filter, @parms) = $sap->ComputeFeatureFilter($source);
    # Determine whether we are looking for protein equivalence or feature
    # equivalence.
    my $protFlag = ($args->{-protein} ? 1 : 0);
    Trace("Protein flag is $protFlag.") if T(3);
    # Loop through the IDs.
    for my $id (@$ids) {
        Trace("Retrieving features for identifier $id.") if T(3);
        # We'll put the FIG IDs for this identifier in here.
        my %fidsFound;
        # Are we looking for proteins or genes?
        if (! $protFlag) {
            # Genes are fairly simple. We use the feature filter to get feature IDs.
            %fidsFound = map { $_ => 1 } $sap->GetFlat($objects, $filter, [@parms, $id],
                                                       'Feature(id)');
        } else {
            # Here we're looking for proteins. We first need to check for a
            # direct protein ID. How we do this depends on the source.
            my ($filter2, @parms2);
            if ($source eq 'mixed') {
                $filter2 = "Identifier(natural-form) = ?";
            } elsif ($source eq 'prefixed') {
                $filter2 = "Identifier(id) = ?"
            } else {
                $filter2 = "Identifier(source) = ? AND Identifier(natural-form) = ?";
                @parms2 = $source;
            }
            my @prots = $sap->GetFlat("Identifier Names ProteinSequence",
                                      $filter2, [@parms2, $id], 'ProteinSequence(id)');
            # Check for proteins related to feature IDs.
            push @prots, $sap->GetFlat("$objects Produces ProteinSequence",
                                       $filter, [@parms, $id], 'ProteinSequence(id)');
            # Now find all the features for these proteins.
            for my $prot (@prots) {
                my @fids = $sap->GetFlat("ProteinSequence IsProteinFor Feature",
                                         "ProteinSequence(id) = ?", [$prot],
                                         'Feature(id)');
                for my $fid (@fids) {
                    $fidsFound{$fid} = 1;
                }
            }
        }
        # Now we apply the genome filter.
        my @results = $sap->FilterByGenome([ keys %fidsFound ], $args->{-genomeName});
        # Put the IDs found into the return hash.
        $retVal->{$id} = [ sort { SeedUtils::by_fig_id($a, $b) } @results ];
    }
    # Return the result.
    return $retVal;
}


=head3 ids_to_genomes

    my $featureHash =       $sapObject->ids_to_genomes({
                                -ids => [$id1, $id2, ...],
                                -source => 'SwissProt',
                                -name => 1
                            });

Return the genome information for each incoming gene ID.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of gene IDs.

=item -source (optional)

Database source of the IDs specified-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
Use C<mixed> to allow mixed ID types (though this may cause problems when the same
ID has different meanings in different databases). Use C<prefixed> to allow IDs with
prefixing indicating the ID type (e.g. C<uni|P00934> for a UniProt ID, C<gi|135813> for
an NCBI identifier, and so forth). The default is C<SEED>.

=item -name (optional)

If TRUE, the genomes names will be returned; if FALSE, the genome IDs will be
returned. The default is FALSE.

=back

=item RETURN

Returns a reference to a hash mapping each incoming ID to the associated genome
ID, or alternatively to the associated genome name.

    $featureHash = { $id1 => $genome1, $id2 => $genome2, ... };

=back

=cut

sub ids_to_genomes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Compute the ID conversion query.
    my ($idObjects, $idFilter, @idParms) = $sap->ComputeFeatureFilter($args->{-source});
    # Determine the desired output field: genome ID or name.
    my $field = ($args->{-name} ? 'Genome(scientific-name)' : 'Genome(id)');
    # Get the incoming IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the incoming IDs.
    for my $id (@$ids) {
        # Get the genome information for this ID.
        my ($genomeID) = $sap->GetFlat("$idObjects IsOwnedBy Genome", $idFilter,
                                       [@idParms, $id], $field);
        # Store it in the return hash.
        $retVal->{$id} = $genomeID;
    }
    # Return the result.
    return $retVal;
}

sub _ids_to_genomes_opt1
{
    my($self, $ids) = @_;

    my $out = $self->_memcache_accelerate($ids, "f", sub {
	my($self, $id_hash, $out, $upd) = @_;

	my @ids = keys %$id_hash;
	my $qs = join(", ", map { "?" } 0..$#ids);
	my $res = $self->{db}->{_dbh}->SQL(qq(SELECT id, function
					      FROM Feature
					      WHERE id IN ($qs)), undef, @ids);
	for my $ent (@$res)
	{
	    my($id, $fn) = @$ent;
	    $out->{$id} = $fn;
	    push(@$upd, ["f:$id", $fn, 12 * 60 * 60]) if $upd;
	}
    });
    return $out;
}

=head3 ids_to_lengths

    my $geneHash =          $sapObjects->ids_to_lengths({
                                -ids => [$id1, $id2, ...],
                                -protein => 1,
                                -source => 'NCBI'
                            });

Return the DNA or protein length of each specified gene.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of gene IDs.

=item -source (optional)

Database source of the IDs specified-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
Use C<mixed> to allow mixed ID types (though this may cause problems when the same
ID has different meanings in different databases). Use C<prefixed> to allow IDs with
prefixing indicating the ID type (e.g. C<uni|P00934> for a UniProt ID, C<gi|135813> for
an NCBI identifier, and so forth). The default is C<SEED>.

=item -genome (optional)

ID of a specific genome. If specified, results will only be returned for genes in the
specified genome. The default is to return results for genes in all genomes.

=item -protein (optional)

If TRUE, then the length of each gene's protein will be returned. Otherwise, the
DNA length of each gene will be returned. The default is FALSE (DNA lengths).

=back

=item RETURN

Returns a reference to a hash mapping each incoming ID to the length of the
associated gene. If no gene is found, or B<-protein> is TRUE and the gene is
not a protein-encoding gene, the ID will not be present in the return hash.

    $geneHash = { $id1 => $length1, $id2 => $length2, ... };

=back

=cut

sub ids_to_lengths {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of feature IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Create the feature filter for IDs of this type.
    my ($objects, $filter, @parms) = $sap->ComputeFeatureFilter($args->{-source},
                                                                $args->{-genome});
    # Determine whether we are looking for DNA or protein lengths.
    my $proteinMode = $args->{-protein} || 0;
    # Loop through the identifiers.
    for my $id (@$ids) {
        # The list of lengths found will be put in here.
        my @lengths;
        # Are we looking for DNA lengths?
        if (! $proteinMode) {
            # Yes. The query is very simple in that case.
            @lengths = $sap->GetFlat($objects, $filter, [@parms, $id], 'Feature(sequence-length)');
        } else {
            # No. We have to get the actual proteins and compute their lengths.
            @lengths = map { length $_ } $sap->GetFlat("$objects Produces ProteinSequence",
                                                       $filter, [@parms, $id],
                                                       'ProteinSequence(sequence)');
        }
        # If we found results, compute the mean length. Most of the time there will only
        # be one result, but some IDs have multiple targets.
        if (@lengths) {
            my $total = 0;
            for my $length (@lengths) { $total += $length; }
            $retVal->{$id} = int($total / @lengths);
        }
    }
    # Return the result hash.
    return $retVal;
}


=head3 make_runs

    my $groupHash =         $sapObject->make_runs({
                                -groups => ["$fid0a, $fid0b, ...",
                                            "$fid1a, $fid1b, ...",
                                            ...],
                                -maxGap => 200,
                                -justFirst = 1,
                                -operonSize => 10000
                            });

Look at sequences of feature IDs and separate them into operons. An
operon contains features that are close together on the same contig going
in the same direction.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -groups

Reference to a list of strings. Each string will contain a comma-separated list of
FIG feature IDs for the features in a group. Alternatively, this can be a
reference to a list of lists, in which each sub-list contains the feature IDs in
a group.

=item -maxGap (optional)

Maximum number of base pairs that can be between to genes in order for them
to be considered as part of the same operon. The default is 200.

=item -justFirst (optional)

If TRUE, then only the first feature in an operon will be included in the
output operon strings. The default is FALSE.

=item -operonSize (optional)

Estimate of the typical size of an operon. This is a tuning parameter; the
default is C<10000>.

=back

=item RETURN

Returns a hash mapping group numbers to lists of operons. In other words,
for each incoming group, the hash will map the group's (zero-based) index number
to a list of operon strings. Each operon string is a comma-separated list of
feature IDs in operon order.

    $groupHash = { 0 => [[$fid1op1, $fid2op1, ...],
                         [$fid1op2, $fid2op2, ...], ... ],
                   1 => [[$fid1opA, $fid2opB, ...],
                         [$fid1opB, $fid2opB, ...], ... ],
                   ... };

=back

=cut

sub make_runs {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of group strings.
    my $groups = ServerThing::GetIdList(-groups => $args);
    my $groupCount = scalar @$groups;
    # Get the just-first flag.
    my $justFirst = $args->{-justFirst} || 0;
    # Get the maximum gap size.
    my $maxGap = $args->{-maxGap} || 200;
    # Get the operon size.
    my $operonSize = $args->{-operonSize} || 10000;
    # Loop through the groups. We use an index because it is going to be
    # the key of the output hash.
    for (my $gidx = 0; $gidx < $groupCount; $gidx++) {
        Trace("Processing group $gidx.") if T(3);
        # Get the features in this group. We allow the option of a comma-delimited
        # string or a sub-list.
        my $group = $groups->[$gidx];
        if (ref $group ne 'ARRAY') {
            $group = [split /\s*,\s*/, $group];
        }
        # This hash is used to remove duplicates.
        my %fids = map { $_ => 1 } @$group;
        # We have our initial set of features. We are now going to loop through
        # all the locations and convert them into operons. Each operon is a list
        # of features. We will collect the completed operons in a list. If at
        # any point we find a feature that's been in a previous operon, we quit.
        # A hash is used to contain the features previously found.
        my %oldFidHash;
        my @operons;
        # This is the main location loop.
        for my $loc (map { $sap->GetLocations($_) } keys %fids) {
            Trace("Computing operon for " . $loc->String . ".") if T(3);
            # We query the database and loop through the locations found
            # until we find a feature from an old operon (which means this
            # operon is discarded) or encounter a gap (which means this operon
            # is complete). We need to go in both directions; the following
            # flag will stop us if it is set in either direction.
            my $redundant = 0;
            # Locations found will be stored in here.
            my @operonData;
            # Search to the left for a gap.
            push @operonData, $sap->FindGapLeft($loc, $maxGap, $operonSize,
                                                \%oldFidHash, \$redundant);
            # Search to the right for a gap.
            push @operonData, $sap->FindGapRight($loc, $maxGap, $operonSize,
                                                 \%oldFidHash, \$redundant);
            # Only proceed if what we found was not redundant.
            if (! $redundant) {
                Trace("Nonredundant operon found.") if T(3);
                # We need to sort the features found into operon order.
                # For a forward operon, we want to sort by leftmost
                # point; for a backward operon, we want to do a reverse
                # sort by rightmost point. The following loop creates
                # a hash we can use for a sort by value.
                my %sortHash;
                for my $operonDatum (@operonData) {
                    # Get the location data from the current tuple.
                    my ($fid, $begin, $dir, $len) = @$operonDatum;
                    # The sort value is something we want to favor when its
                    # numerically low. For a forward operon, this is the
                    # begin point. For a backward operon, we take the end point
                    # and negate it. The negation makes it sort correctly.
                    my $sortValue = ($dir eq '+' ? $begin : -($begin + $len));
                    # Merge this datum into the hash.
                    if (! exists $sortHash{$fid}) {
                        # Here we have a new feature. Save its sort value.
                        $sortHash{$fid} = $sortValue;
                        # Add it to the redundancy hash for future use.
                        $oldFidHash{$fid} = 1;
                    } else {
                        # Hewre we have a second location for an existing feature.
                        # Save the minimum sort value.
                        $sortHash{$fid} = Tracer::Min($sortValue, $sortHash{$fid});
                    }
                }
                # Create the operon.
                my @operon = sort { $sortHash{$a} <=> $sortHash{$b} } keys %sortHash;
                # Add it to the return list, in the format indicated by the
                # "justFirst" flag.
                if ($justFirst) {
                    push @operons, $operon[0];
                } else {
                    push @operons, join(", ", @operon);
                }
            }
        }
        # Put the operons found into the return hash.
        $retVal->{$gidx} = \@operons;
    }
    # Return the result.
    return $retVal;
}

=head3 proteins_to_fids

    my $protHash =          $sapObject->proteins_to_fids({
                                -prots => [$prot1, $prot2, ...]
                            });

Return the FIG feature IDs associated with each incoming protein. The protein can be
specified as an amino acid sequence or MD5 protein ID.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -prots

Reference to a list of proteins. Each protein can be specified as either an amino acid
sequence or an MD5 protein ID. The method will assume a sequence of 32 hex characters is
an MD5 ID and anything else is an amino acid sequence. Amino acid sequences should be
in upper-case only.

=back

=item RETURN

Returns a hash mapping each incoming protein to a list of FIG feature IDs for the genes that
produce the protein.

    $protHash = { $prot1 => [$fid1a, $fid1b, ...],
                  $prot2 => [$fid2a, $fid2b, ...],
                  ... };

=back

=cut

sub proteins_to_fids {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of incoming proteins.
    my $prots = ServerThing::GetIdList(-prots => $args);
    # Declare the return hash.
    my $retVal = {};
    # Loop through the proteins.
    my @ids;
    my $have_prot;
    for my $prot (@$prots) {
        # If this is a protein sequence, convert it to an MD5 ID.
        my $id;
        if (length($prot) == 32 && $prot !~ /[^0-9a-f]/) {
            $id = $prot;
        } else {
            $id = $sap->ProteinID($prot);
	    $have_prot++;
        }
	push(@ids, [$prot, $id]);
    }

    if ($self->{memcache} && !$have_prot)
    {
	my $opt_out = $self->_proteins_to_fids_opt1($prots);
	if ($opt_out->{$prots->[0]} ne '')
	{
	    return $opt_out;
	}
	else
	{
	    print STDERR "proteins_to_fids: missing output " . Dumper($prots, $opt_out);
	}
    }

    for my $id_ent (@ids)
    {
	my($prot, $id) = @$id_ent;
        # Find the features that produce the specified protein.
        my @fids = $sap->GetFlat('IsProteinFor', 'IsProteinFor(from-link) = ?', [$id], "to-link");
        # Put the results in the return hash.
        $retVal->{$prot} = \@fids;
    }
    # Return the results.
    return $retVal;
}

sub _proteins_to_fids_opt1
{
    my($self, $ids) = @_;

    my $out = $self->_memcache_accelerate_list($ids, "md52f", sub {
	my($self, $id_hash, $out, $upd) = @_;

	my @ids = keys %$id_hash;
	my $qs = join(", ", map { "?" } 0..$#ids);
	my $res = $self->{db}->{_dbh}->SQL(qq(SELECT from_link, to_link
					      FROM IsProteinFor
					      WHERE from_link IN ($qs)), undef, @ids);
	for my $ent (@$res)
	{
	    my($md5, $id) = @$ent;
	    push(@{$out->{$md5}}, $id);
	    push(@{$upd->{"md52f:$md5"}}, $id);
	}

    });
    return $out;
}

=head2 FIGfam Data Methods

=head3 all_figfams

    my $ffHash =            $sapObject->all_figfams({
                                -roles => [$role1, $role2, ...],
                                -functions => [$function1, $function2, ...]
    });

Return a list of all the FIGfams along with their functions. Optionally, you
can specify a role or a function, and only FIGfams with that role or function
will be returned.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item roles (optional)

If specified, a reference to a list of roles. Only FIGfams with one of the
specified roles (or one of the functions listed in C<-functions>) will be
returned in the hash.

=item function (optional)

If specified, a reference to a list of functions. Only FIGfams with one of the
specified functions (or one of the roles listed in C<-roles>) will be returned
in the hash.

=back

=item RETURN

Returns a reference to a hash mapping each qualifying FIGfam ID to its
function.

    $ffHash = { $ff1 => $function1, $ff2 => $function2, ... };

=back

=cut

sub all_figfams {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the roles and functions. Both of these are optional.
    my $roles = ServerThing::GetIdList(-roles => $args, 1);
    my $functions = ServerThing::GetIdList(-functions => $args, 1);
    # If we have neither, we're officially asking for everything.
    if (@$roles + @$functions == 0) {
        # Ask for all of the FIGfams.
        $retVal = { map { $_->[0] => $_->[1] } $sap->GetAll('Family',
                                                            'Family(id) LIKE ?',
                                                            ['FIG%'],
                                                            ['id', 'family-function']) };
    } else {
        # Here we are searching by role or function. Create hashes for the
        # roles and functions of interest.
        my %roleMap = map { $_ => 1 } @$roles;
        my %functionMap = map { $_ => 1 } @$functions;
        # Now we want a list of all roles of interest. We start with the roles
        # themselves, then add all the roles taken from functions.
        my %roleList = map { $_ => 1 } @$roles;
        for my $function (@$functions) {
            for my $role (roles_of_function($function)) {
                $roleList{$role} = 1;
            }
        }
        # Now loop through all of the roles found.
        for my $role (keys %roleList) {
            # Get all the FIGfams for this role.
            my @ffPairs = $sap->GetAll("Role DeterminesFunctionOf Family",
                                       'Role(id) = ? AND Family(id) LIKE ?', [$role,'FIG%'],
                                       [qw(Family(id) Family(family-function))]);
            # Loop through the FIGfams for this role.
            for my $ffPair (@ffPairs) {
                my ($ff, $function) = @$ffPair;
                # Keep this FIGfam if the role or function is of interest.
                if ($roleMap{$role} || $functionMap{$function}) {
                    $retVal->{$ff} = $function;
                }
            }
        }
    }
    # Return the result.
    return $retVal;
}

=head3 discriminating_figfams

    my $groupList =         $sapObject->discriminating_figfams({
                                -group1 => [$genome1a, $genome2a, ...],
                                -group2 => [$genome2a, $genome2b, ...]
                            });

Determine the FIGfams that discriminate between two groups of genomes.

A FIGfam discriminates between genome groups if it is common in one group and
uncommon in the other. The degree of discrimination is assigned a score based
on statistical significance, with 0 being insignificant and 2 being extremely
significant. FIGfams with a score greater than 1 are returned by this method.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -group1

Reference to a list of genome IDs for the genomes in the first group.

=item -group2

Reference to a list of genome IDs for the genomes in the second

=back

=item RETURN

Returns a reference to a 2-tuple, consisting of (0) a hash mapping FIGfam IDs
to scores for FIGfams common in group 1 and (1) a hash maping FIGfam IDs to
scores for FIGfams common in group 2.

    $groupList = [{ $ff1a => $score1a, $ff1b => $score1b, ... },
                  { $ff2a => $score2a, $ff2b => $score2b, ... }];

=back

=cut

sub discriminating_figfams {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the two groups.
    my $group1 = ServerThing::GetIdList(-group1 => $args);
    my $group2 = ServerThing::GetIdList(-group2 => $args);
    # Generate the FIGfam hashes for the groups.
    my $group1H = $self->genome_figfams({ -ids => $group1 });
    my $group2H = $self->genome_figfams({ -ids => $group2 });
    # Get the signature tool.
    require Signatures;
    # Compute the signature.
    my ($sig1H, $sig2H) = Signatures::ComputeSignatures($group1H, $group2H);
    # Return the result.
    return [$sig1H, $sig2H];
}

=head3 figfam_fids

    my $fidList =           $sapObject->figfam_fids({
                                -id => $figFam1,
                                -fasta => 1
                            });

Return a list of all the protein encoding genes in a FIGfam. The genes
can be returned as IDs or as FASTA strings.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -id

ID of the desired FIGfam.

=item -fasta

If TRUE, then the output will be in the form of FASTA strings; otherwise it will
be in the form of FIG IDs.

=back

=item RETURN

Returns a reference to a list of genes in the form of FIG feature IDs or protein
FASTA strings.

=over 8

=item Normal Output

    $fidList = [$fid1, $fid2, ...];

=item Output When -fasta = 1

    $fidList = [$fasta1, $fasta2, ...];

=back

=back

=cut

sub figfam_fids {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = [];
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the incoming FIGfam ID.
    my $figFam = $args->{-id};
    Confess("No FIGfam ID specified.") if ! defined $figFam;
    # Get the FASTA-format flag.
    my $fasta = $args->{-fasta} || 0;
    # Configure the query depending on whether or not we're doing the FASTA
    # sequences.
    my ($objects, $fields) = ("HasMember", ["HasMember(to-link)"]);
    if ($fasta) {
        $objects .= " Feature Produces ProteinSequence";
        push @$fields, "HasMember(from-link)", "ProteinSequence(sequence)";
    }
    # Get the list of genes in the FIGfam.
    my @fids = $sap->GetAll($objects, "HasMember(from-link) = ?", [$figFam],
                            $fields);
    # Are we doing FASTA results?
    if (! $fasta) {
        # No, so just return the FIDs.
        $retVal = [ map { $_->[0] } @fids ];
    } else {
        # Yes, so we need a FASTA string for each feature.
        $retVal = [ map { create_fasta_record(@$_) } @fids ];
    }
    # Return the result.
    return $retVal;
}

=head3 figfam_fids_batch

    my $fidHash =           $sapObject->figfam_fids_batch({
                                -ids => [$ff1, $ff2, ...],
                                -genomeFilter => $genome1
                            });

Return a list of all the protein encoding genes in one or more FIGfams. This
method is an alternative to L</figfam_fids> that is faster when you need the
feature IDs but not the protein sequences.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the IDs of the desired FIGfams.

=item -genomeFilter (optional)

The ID of a genome. If specified, then only feature IDs from the specified
genome will be returned.

=back

=item RETURN

Returns a hash mapping each incoming FIGfam ID to a list of the IDs for the features
in that FIGfam.

    $fidHash = { $ff1 => [$fid1a, $fid1b, ...],
                 $ff2 => [$fid2a, $fid2b, ...],
                 ... };

=back

=cut

sub figfam_fids_batch {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of FIGfam IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    my %figFamHash = map { $_ => [] } @$ids;
    Trace("Initiating batch FIGfam retrieval for " . scalar(@$ids) . " FIGfams.") if T(3);
    # Form the filter string. At a minimum, we filter by FIGfam ID.
    my $filter = "HasMember(from-link) = ?";
    my @parms;
    # If there is a genome filter, incorporate that.
    my $genomeFilter = $args->{-genomeFilter};
    if ($genomeFilter) {
        $filter .= " AND HasMember(to-link) LIKE ?";
        push @parms, "fig|$genomeFilter.%";
    }
    # Loop through the FIGfams.
    for my $id (sort keys %figFamHash) {
        # Get the list of genes in the FIGfam.
        my @fids = $sap->GetFlat('HasMember', $filter, [$id, @parms],
                                'to-link');
        # Store them in the hash.
        $figFamHash{$id} = \@fids;
    }
    Trace("Batch FIGfam retrieval complete.") if T(3);
    # Return the result.
    return \%figFamHash;
}

=head3 figfam_function

    my $ffHash =            $sapObject->figfam_function({
                                -ids => [$ff1, $ff2, ...]
                            });

For each incoming FIGfam ID, return its function, that is, the common
functional assignment of all its members.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIGfam IDs.

=back

=item RETURN

Returns a hash mapping each incoming FIGfam ID its function string.

    $ffHash => { $ff1 => $function1, $ff2 => $function2, ... };

=back

=cut

sub figfam_function {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of FIGfam IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the FIGfams.
    for my $id (@$ids) {
        # Get this FIGfam's function.
        my ($function) = $sap->GetEntityValues(Family => $id, ['family-function']);
        # Put it in the return hash.
        $retVal->{$id} = $function;
    }
    # Return the result.
    return $retVal;
}

=head3 genome_figfams

    my $genomeHash =        $sapObject->genome_figfams({
                                -ids => [$genome1, $genome2, ...]
                            });

Compute the list of FIGfams represented in each specific genome.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of genome identifiers.

=back

=item RETURN

Returns a reference to a hash mapping each incoming genome ID to a list of the
IDs of the FIGfams represented in that genome.

    $genomeHash = { $genome1 => [$ff1a, $ff1b, ...],
                    $genome2 => [$ff2a, $ff2b, ...],
                     ... };

=back

=cut

sub genome_figfams {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return hash.
    my $retVal = {};
    # Get the list of genomes.
    my $genomes = ServerThing::GetIdList(-ids => $args);
    # Loop through the genomes, finding FIGfams.
    for my $genome (@$genomes) {
        # Get this genome's FIGfam list.
        my @figfams = $sap->GetFlat("HasRepresentativeOf",
                                     'HasRepresentativeOf(from-link) = ?', [$genome],
                                     'to-link');
        # Store it in the return hash.
        $retVal->{$genome} = \@figfams;
    }
    # Return the result.
    return $retVal;
}

=head3 ids_to_figfams

    my $featureHash =       $sapObject->ids_to_figfams({
                                -ids => [$id1, $id2, ...],
                                -functions => 1,
                                -source => 'RefSeq'
                            });

This method returns a hash mapping each incoming feature to its FIGfam.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of feature identifiers.

=item -functions (optional)

If TRUE, the family function will be returned in addition to the list of
FIGfam IDs. In this case, instead of a list of FIGfam IDs, each feature ID will
point to a list of 2-tuples, each consisting of (0) a FIGfam ID followed by (1)
a function string. The default is FALSE.

=item -source (optional)

Database source of the IDs specified-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those
databases. Use C<mixed> to allow mixed ID types (though this may cause problems
when the same ID has different meanings in different databases). Use C<prefixed>
to allow IDs with prefixing indicating the ID type (e.g. C<uni|P00934> for a
UniProt ID, C<gi|135813> for an NCBI identifier, and so forth). The default is
C<SEED>.

=item -genome (optional)

ID of a specific genome. If specified, results will only be returned for genes in the
specified genome. The default is to return results for genes in all genomes.

=back

=item RETURN

Returns a reference to a hash mapping each incoming feature ID to a list of the
IDs of the FIGfams that contain it. (In general the list will be a singleton
unless the feature ID corresponds to multiple actual features.) Features not in
FIGfams will be omitted from the hash.

    $featureHash = { $id1 => [$ff1a, $ff1b, ...],
                     $id2 => [$ff2a, $ff2b, ...],
                     ... };

=back

=cut

sub ids_to_figfams {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the feature ID list.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Find out if we need to see functions in the output.
    my $functionFlag = $args->{-functions} || 0;
    # Compute the ID conversion query.
    my ($idObjects, $idFilter, @idParms) = $sap->ComputeFeatureFilter($args->{-source},
                                                                      $args->{-genome});
    # Loop through the incoming feature IDs, retrieving FIGfams.
    for my $id (@$ids) {
        Trace("Reading families for $id.") if T(3);
        # Get this feature's FIG families and functions.
        my @fidPairs = $sap->GetAll("$idObjects IsMemberOf Family",
                                    "$idFilter AND Family(id) LIKE ?",
                                    [@idParms, $id, 'FIG%'],
                                    ['Family(id)', 'Family(family-function)']);
        # Only proceed if we found something.
        if (@fidPairs) {
            # Do we want functions or just FIGfams?
            if ($functionFlag) {
                # We want both.
                $retVal->{$id} = \@fidPairs;
            } else {
                # We want only the FIGfams.
                $retVal->{$id} = [ map { $_->[0] } @fidPairs ];
            }
        }
    }
    # Return the result.
    return $retVal;
}

=head3 related_figfams

    my $ffHash =            $sapObject->related_figfams({
                                -ids => [$ff1, $ff2, ...],
                                -expscore => 1,
                                -all => 1
                            });

This method takes a list of FIGfam IDs. For each FIGfam, it returns a
list of FIGfams related to it by functional coupling.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIGfam IDs.

=item -expscore (optional)

If TRUE, then the score returned will be the co-expression score. If
FALSE, the score returned will be the co-occurrence score. This option
is ignored if C<-all> is specified. The default is FALSE.

=item -all (optional)

If TRUE, then both scores will be returned. The default is FALSE, meaning
only one score is returned.

=back

=item RETURN

=over 8

=item normal

Returns a reference to a hash mapping each incoming FIGfam ID
to a list of 2-tuples for other FIGfams. The 2-tuples
each consist of (0) a related FIGfam's ID followed by (1) a 2-tuple
containing a coupling score and the related FIGfam's function.

    $ffHash = { $ff1 => [[$ff1a, [$score1a, $function1a]],
                         [$ff1b, [$score1b, $function1b]], ...],
                $ff2 => [[$ff2a, [$score2a, $function2a]],
                         [$ff2b, [$score2b, $function2b]], ...],
                ... };

=item -exp = all

Returns a reference to a hash mapping each incoming FIGfam ID
to a list of 2-tuples for other FIGfams. The 2-tuples
each consist of (0) a related FIGfam's ID followed by (1) a 3-tuple
containing the co-occurrence coupling score, the co-expression coupling
score, and the related FIGfam's function.

    $ffHash = { $ff1 => [[$ff1a, [$score1ax, $score1ay, $function1a]],
                         [$ff1b, [$score1bx, $score1by, $function1b]], ...],
                $ff2 => [[$ff2a, [$score2ax, $score2ay, $function2a]],
                         [$ff2b, [$score2bx, $score2by, $function2b]], ...],
                ... };


=back

=back

=cut

sub related_figfams {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of FIGfam IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Determine the score that we want.
    my @scoreTypes;
    if ($args->{-all}) {
        @scoreTypes = qw(co-occurrence-evidence co-expression-evidence);
    } elsif ($args->{-expscore}) {
        @scoreTypes = qw(co-expression-evidence);
    } else {
        @scoreTypes = qw(co-occurrence-evidence);
    }
    # Loop through the FIGfams.
    for my $id (@$ids) {
        # The couplings for this family will be put in here.
        my @couples;
        # We need to do this query for both directions of the coupling relationship.
        for my $rel (qw(IsCoupledWith IsCoupledTo)) {
            # Get a query to find the familes coupled in this direction.
            my $qh = $sap->Get("$rel Family", "$rel(from-link) = ?", [$id]);
            while (my $resultRow = $qh->Fetch()) {
                # Get this coupled family.
                my $familyId = $resultRow->PrimaryValue('Family(id)');
                my $familyFamilyFunction = $resultRow->PrimaryValue('Family(family-function)');
                # Get the scores.
                my @scores;
                for my $scoreType (@scoreTypes) {
                    push @scores, $resultRow->PrimaryValue("$rel($scoreType)");
                }
                # Store this coupling in the current output list.
                push @couples, [$familyId, [@scores, $familyFamilyFunction]];
            }
        }
        # Store the list of results.
        $retVal->{$id} = \@couples
    }
    # Return the result.
    return $retVal;
}

=head3 roles_to_figfams

    my $roleHash =          $sapObject->roles_to_figfams({
                                -roles => [$role1, $role2, ...]
                            });

For each incoming role, return a list of the FIGfams that implement
the role, that is, whose functional assignments include the role.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -roles

Reference to a list of role names.

=back

=item RETURN

Returns a reference to a hash mapping each incoming role to a list of
FIGfam IDs for the FIGfams that implement the role.

    $roleHash = { $role1 => [$ff1a, $ff1b, ...],
                  $role2 => [$ff2a, $ff2b, ...],
                  ... };

=back

=cut

sub roles_to_figfams {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the list of roles.
    my $roles = ServerThing::GetIdList(-roles => $args);
    # Declare the return hash.
    my $retVal = {};
    # Loop through the roles.
    for my $role (@$roles) {
        # Get the FIGfams for this role.
        my @ffs = $sap->GetFlat("DeterminesFunctionOf",
            'DeterminesFunctionOf(from-link) = ?', [$role],
            'to-link');
        # Store them in the return hash.
        $retVal->{$role} = \@ffs;
    }
    # Return the result hash.
    return $retVal;
}

=head2 Functional Coupling Data Methods

=head3 clusters_containing

    my $featureHash =       $sapObject->clusters_containing({
                                -ids => [$fid1, $fid2, ...]
                            });

This method takes as input a list of FIG feature IDs. For each feature, it
returns the IDs and functions of other features in the same cluster of
functionally-coupled features.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs.

=back

For backward compatibility, this method can also take as input a reference to
a list of FIG feature IDs.

=item RETURN

Returns a reference to a hash. The hash maps each incoming feature ID to a
2-tuple containing (0) the feature's functional assignment and (1) a
reference to a hash that maps each clustered feature to its functional assignment.

    $featureHash = { $fid1 => [$function1, { $fid1a => $function1a,
                                             $fid1b => $function1b,
                                             ...}],
                     $fid2 => [$function2, { $fid2a => $function2a,
                                             $fid2b => $function2b,
                                             ...}],
                     ... };

In backward-compatibility mode, this method returns a reference to a list. For
each incoming feature, there is a list entry containing the feature ID, the
feature's functional assignment, and a sub-list of 2-tuples. Each 2-tuple
contains the ID of another feature in the same cluster and its functional
assignment.

=back

=cut

sub clusters_containing {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get access to functional coupling services.
    require FC;
    # Check for backward-compatibility mode.
    my $backwardMode = 0;
    if (ref $args ne 'HASH') {
        $args = { -ids => $args };
        $backwardMode = 1;
    }
    # Get the list of feature IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the features.
    for my $id (@$ids) {
        # Get this feature's cluster data.
        my $cluster = &FC::in_co_occurrence_cluster($sapling, $id);
        # Did we find something?
        if ($cluster) {
            # Get this feature's assignment.
            my $func = scalar $sapling->Assignment($id);
            # Create a hash of the clustered IDs.
            my %members = map { $_ => $sapling->Assignment($_) } @$cluster;
            # Store the result.
            $retVal->{$id} = [$func, \%members];
        }
    }
    # In backward-compatibility mode, convert the result to a list.
    if ($backwardMode) {
        # We'll create our result list in here.
        my @outList;
        # Loop through the IDs.
        for my $id (@$ids) {
            # Do we have something for this feature?
            my $featureData = $retVal->{$id};
            if (defined $featureData) {
                # Get the pieces.
                my ($func, $memberHash) = @$featureData;
                # Convert the member hash to a list of 2-tuples.
                my @memberList = map { [$_, $memberHash->{$_} ] } sort keys %$memberHash;
                # Assemble the result.
                push @outList, [$id, $func, \@memberList];
            }
        }
        # Store the reformatted list.
        $retVal = \@outList;
    }
    # Return the result.
    return $retVal;
}

=head3 co_occurrence_evidence

    my $pairHash =          $sapObject->co_occurrence_evidence({
                                -pairs => ["$fid1:$fid2", "$fid3:$fid4", ...]
                            });

For each specified pair of genes, this method returns the evidence that
the genes are functionally coupled (if any); that is, it returns a list
of the physically close homologs for the pair.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -pairs

Reference to a list of functionally-coupled pairs. Each pair is represented by two
FIG gene IDs, either in the form of a 2-tuple or as a string with the two gene IDs
separated by a colon.

=back

=item RETURN

Returns a hash mapping each incoming gene pair to a list of 2-tuples. Each 2-tuple
contains a pair of physically close genes, the first of which is similar to the first
gene in the input pair, and the second of which is similar to the second gene in the
input pair. The hash keys will consist of the two gene IDs separated by a colon (e.g.
C<fig|273035.4.peg.1016:fig|273035.4.peg.1018>).

    $pairHash = { "$fid1:$fid2" => [[$fid1a, $fid2a], [$fid1b, $fid2b], ...],
                  "$fid3:$fid4" => [[$fid3a, $fid4a], [$fid3b, $fid4b], ...],
                  ... };

=back

=cut

sub co_occurrence_evidence {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get access to the functional coupling services.
    require FC;
    # Get the list of pairs.
    my $pairs = ServerThing::GetIdList(-pairs => $args);
    # Loop through the pairs.
    for my $pair (@$pairs) {
        # Determine the IDs in this pair.
        my ($peg1, $peg2);
        if (ref $pair) {
            ($peg1, $peg2) = @$pair;
        } else {
            ($peg1, $peg2) = split /:/, $pair;
        }
        # Get the evidence and store it in the return hash.
        $retVal->{"$peg1:$peg2"} = FC::co_occurrence_evidence($sap, $peg1, $peg2);
    }
    # Return the result.
    return $retVal;
}

=head3 conserved_in_neighborhood

    my $featureHash =       $sapObject->conserved_in_neighborhood({
                                -ids => [$fid1, $fid2, ...]
                            });

This method takes a list of feature IDs. For each feature ID, it will
return the set of other features to which it is functionally coupled,
along with the appropriate score.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs.

=back

For backward compatibility, this method can also take as input a reference to
a list of FIG feature IDs.

=item RETURN

Returns a reference to a hash mapping each incoming feature ID to a list of
4-tuples, one 4-tuple for each feature coupled to the incoming feature. Each
4-tuple contains (0) the coupling score, (1) the FIG ID of the coupled feature,
(2) the coupled feature's current functional assignment, and (3) the ID of the
pair set to which the coupling belongs.

    $featureHash = { $fid1 => [[$score1A, $fid1A, $function1A, $psID1A],
                               [$score1B, $fid1B, $function1B, $psID1B], ...],
                     $fid2 => [[$score2A, $fid2A, $function2A, $psID2A],
                               [$score2B, $fid2B, $function2B, $psID2B], ...],
                     ... };

In backward compatibility mode, returns a list of sub-lists, each sub-list
corresponding to the value that would be found in the hash for the feature in the
specified position of the input list.

=back

=cut

sub conserved_in_neighborhood {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # Get the functional coupling methods.
    require FC;
    # Declare the return variable.
    my $retVal = {};
    # Check for backward compatibility mode.
    my $backwardMode = 0;
    # Convert a list to a hash.
    if (ref $args ne 'HASH') {
        $args = { -ids => $args };
        $backwardMode = 1;
    }
    # Get the list of feature IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the features.
    for my $id (@$ids) {
        # Create a sub-list for this feature.
        my $group = [];
        # Ask for the functional coupling information.
        my @co_occurs = &FC::co_occurs($sapling, $id);
        # Loop through the coupling data found.
        for my $tuple (@co_occurs) {
            # Get the coupled feature's data.
            my($sc, $fid, $pairset) = @$tuple;
            # Add it to the group of tuples for this feature's couplings.
            push(@$group, [$sc, $fid, $sapling->Assignment($fid), $pairset]);
        }
        # Add this feature's couplings to the return value.
        $retVal->{$id} = $group;
    }
    # If we're in backward-compatibility mode, convert the output to a list.
    if ($backwardMode) {
        my @outList = map { $retVal->{$_} } @$ids;
        $retVal = \@outList;
    }
    # Return the result.
    return $retVal;
}

=head3 pairsets

    my $psHash =            $sapObject->pairsets({
                                -ids => [$psID1, $psID2, ...]
                            });

This method takes as input a list of functional-coupling pair set IDs
(such as those returned in the output of L</conserved_in_neighborhood>). For
each pair set, it returns the set's score (number of significant couplings) and
a list of the coupled pairs in the set.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of functional-coupling pair set IDs.

=back

For backward compatibility, you may also specify a reference to a list of pair
set IDs.

=item RETURN

Returns a reference to a hash that maps each incoming pair-set ID to a 2-tuple
that consists of (0) the set's score and (1) a reference to a list of 2-tuples
containing the pairs in the set.

    $psHash = { $psID1 => [$score1, [[$fid1A, $fid1B],
                                     [$fid1C, $fid1D], ...]],
                $psID2 => [$score2, [[$fid2A, $fid2B],
                                     [$fid2C, $fid2D], ...]],
                ... };

In backward-compatibility mode, returns a reference to a list of 2-tuples, each
consisting of (0) an incoming pair-set ID, and (1) the 2-tuple that would be its
hash value in the normal output.

=back

=cut

sub pairsets {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # Get access to the functional coupling methods.
    require FC;
    # Declare the return variable.
    my $retVal = {};
    # Check for backward-compatability mode.
    my $backwardMode = 0;
    if (ref $args ne 'HASH') {
        $args = { -ids => $args };
        $backwardMode = 1;
    }
    # Get the list of pairset IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the pairsets, producing output.
    for my $id (@$ids) {
        $retVal->{$id} = [&FC::co_occurrence_set($sapling, $id)];
    }
    # In backward-compatible mode, convert the output to a list.
    if ($backwardMode) {
        my @outList = map { [$_, $retVal->{$_}] } @$ids;
        $retVal = \@outList;
    }
    # Return the result.
    return $retVal;
}


=head3 related_clusters

    my $featureHash =       $sapObject->related_clusters({
                                -ids => [$fid1, $fid2, ...]
                            });

This method returns the functional-coupling clusters related to the specified
input features. Each cluster contains features on a single genome that are
related by functional coupling.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of FIG feature IDs.

=back

=item RETURN

Returns a reference to a hash that maps each incoming feature ID to a list
of clusters. Each cluster in the list is a 3-tuple consisting of (0) the ID of a
feature similar to the incoming feature, (1) the similarity P-score, and (2) a
reference to a list of 2-tuples containing clustered features and their functional
assignments.

    $featureHash = { $fid1 => [[$fid1A, $score1A, [[$fid1Ax, $function1Ax],
                                                   [$fid1Ay, $function1Ay],
                                                   ...]],
                               [$fid1B, $score1B, [[$fid1Bx, $function1Bx],
                                                   [$fid1By, $function1By],
                                                   ...]],
                               ...],
                      $fid2 => [[$fid2A, $score2A, [[$fid2Ax, $function2Ax],
                                                   [$fid2Ay, $function2Ay],
                                                   ...]],
                               [$fid2B, $score2B, [[$fid2Bx, $function2Bx],
                                                   [$fid2By, $function2By],
                                                   ...]],
                               ...],
                      ... };

=back

=cut

sub related_clusters {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the functional coupling methods.
    require FC;
    # Get the list of feature IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the features.
    for my $id (@$ids) {
        # Create the output list for this feature.
        my $output = [];
        # Loop through the related clusters.
        for my $cluster (FC::largest_co_occurrence_clusters($sapling, $id)) {
            # Get this cluster's data.
            my ($fid, $sc, $other_fids) = @$cluster;
            # Extract the functional roles of the other features in the cluster.
            my $other_tuples = [ map { [$_, $sapling->Assignment($_)] } @$other_fids ];
            # Assemble the result into the output list.
            push @$output, [$fid, $sc, $other_tuples];
        }
        # Return this list of clusters.
        $retVal->{$id} = $output;
    }
    # Return the result.
    return $retVal;
}


=head2 Genome Data Methods

=head3 all_features

    my $genomeHash =        $sapObject->all_features({
                                -ids => [$genome1, $genome2, ...],
                                -type => [$type1, $type2, ...],
                            });

Return a list of the IDs for all features of a specified type in a specified
genome.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of genome IDs.

=item -type (optional)

Type of feature desired (e.g. C<peg>, C<rna>), or a reference to a list of
desired feature types. If omitted, all features regardless of type are returned.

=back

=item RETURN

Returns a reference to a hash that maps each incoming genome ID to a list
of the desired feature IDs for that genome. If a genome does not exist or has
no features of the desired type, its ID will map to an empty list.

    $genomeHash = { $genome1 => [$fid1a, $fid1b, ...],
                    $genome2 => [$fid2a, $fid2b, ...],
                    ... };

=back

=cut

sub all_features {
    # Get the parameters.
    my ($self, $args) = @_;
    # Create the filter for the query.
    my $filter = 'IsOwnerOf(from-link) = ?';
    # The query may require additional parameters in addition to the genome ID.
    # Those additional parameters will go in here.
    my @parms;
    # Are we filtering by type?
    my $type = $args->{-type};
    if (defined $type) {
        # Yes. Add filtering by type. First, insure that we are dealing with a
        # list of types.
        if (ref $type ne 'ARRAY') {
            $type = [$type];
        }
        # Now, form the list into an IN-type filter.
        $filter .= " AND Feature(feature-type) IN (" . join(", ", map { "?" } @$type) . ")";
        push @parms, @$type;
    }
    # Declare the return hash.
    my %retVal;
    # Get the list of genome IDs.
    my $genomeIDs = ServerThing::GetIdList(-ids => $args);
    # Loop through the genome IDs.
    for my $genomeID (@$genomeIDs) {
        # Execute the query.
        Trace("Retrieving features for $genomeID.") if T(3);
        my @fids = $self->{db}->GetFlat("IsOwnerOf Feature", $filter,
                                        [$genomeID, @parms], 'Feature(id)');
        # Put the result in the output hash.
        $retVal{$genomeID} = \@fids;
    }
    # Return the hash of results.
    return \%retVal;
}

=head3 all_genomes

    my $genomeHash = $sapObject->all_genomes({
                            -complete => 1,
                            -prokaryotic => 1
                        });

Return a list of the IDs for all the genomes in the system.

=over 4

=item parameter

Reference to a hash containing the following keys.

=over 8

=item -complete (optional)

If TRUE, only complete genomes will be returned. The default is FALSE (return
all genomes).

=item -prokaryotic (optional)

If TRUE, only prokaryotic genomes will be returned. The default is FALSE (return
all genomes).

=back

=item RETURN

Returns a reference to a hash mapping genome IDs to genome names.

    $genomeHash = { $genome1 => $name1, $genome2 => $name2, ... };

=back

=cut

sub all_genomes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Fix the filter and parms according to the options.
    my @filters;
    my $parms = [];
    if ($args->{-complete}) {
        push @filters, "Genome(complete) = ?";
        push @$parms, 1;
    }
    if ($args->{-prokaryotic}) {
        push @filters, "Genome(prokaryotic) = ?";
        push @$parms, 1;
    }
    my $filter = join(" AND ", @filters);
    # Ask for the genome data
    my %retVal = map { $_->[0] => $_->[1] }
                    $self->{db}->GetAll("Genome", $filter, $parms,
                                        "Genome(id) Genome(scientific-name)");
    # Return the result.
    return \%retVal;
}

=head3 all_proteins

    my $fidHash = $sapObject->all_proteins({
                        -id => $genome1
                    });

Return the protein sequences for all protein-encoding genes in the specified
genome.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -id

A single genome ID. All of the protein sequences for genes in the specified
genome will be extracted.

=back

=item RETURN

Returns a reference to a hash that maps the FIG ID of each protein-encoding
gene in the specified genome to its protein sequence.

    $fidHash = { $fid1 => $protein1, $fid2 => $protein2, ... };

=back

=cut

sub all_proteins {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the target genome ID.
    my $genome = $args->{-id};
    Confess("No genome ID specified for all_proteins.") if ! defined $genome;
    # Ask for all the proteins in this genome and put them in a hash.
    my %retVal = map { $_->[0] => $_->[1] }
                    $sap->GetAll("Produces ProteinSequence",
                                 'Produces(from-link) LIKE ?', ["fig|$genome.peg.%"],
                                 [qw(Produces(from-link) ProteinSequence(sequence))]);
    # Return the result.
    return \%retVal;
}

=head3 close_genomes

    my $genomeHash = $sapObject->close_genomes({
                        -ids => [$genome1, $genome2, ...],
                        -count => 10,
                    });

Find the genomes functionally close to the input genomes.

Functional closeness is determined by the number of FIGfams in common. As a result,
this method will not produce good results for genomes that do not have good FIGfam
coverage.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item ids

Reference to a list of genome IDs for the genomes whose close neighbors are
desired.

=item count (optional)

Maximum number of close genomes to return for each input genome. The default is
C<10>.

=back

=item RETURN

Returns a reference to a hash mapping each incoming genome ID to a list of
2-tuples. Each 2-tuple consists of (0) the ID of a close genome and (2) the
score (from 0 to 1) for the match. The list will be sorted from closest to
furthest.

=back

=cut

sub close_genomes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Compute the count.
    my $count = $args->{-count} || 10;
    # Loop through the IDs, locating the close genomes.
    for my $id (@$ids) {
        # This hash will count the common features in related genomes.
        my %genomes;
        # This hash will count the FIGfams in the source genome.
        my %figFams;
        # Get the figfam data for this genome.
        my $qh = $sap->Get("IsOwnerOf IsMemberOf HasMember",
                           'IsOwnerOf(from-link) = ?', [$id]);
        while (my $resultRow = $qh->Fetch()) {
            # Get this figfam and the related feature ID.
            my $figFam = $resultRow->PrimaryValue('HasMember(from-link)');
            my $fid = $resultRow->PrimaryValue('HasMember(to-link)');
            # Compute the feature's genome.
            my $genome = genome_of($fid);
            # Count this genome.
            $genomes{$genome}++;
            # Count this FIGfam;
            $figFams{$figFam} = 1;
        }
        # Sort the genomes and compute the scores.
        my @genomes = sort { -($genomes{$a} <=> $genomes{$b}) } keys %genomes;
        my $figFamCount = scalar keys %figFams;
        # We'll put our results in here.
        my @results;
        for my $genome (@genomes) { last if scalar(@results) >= $count;
            if ($genome ne $id) {
                push @results, [$genome, $genomes{$genome} / $figFamCount ];
            }
        }
        $retVal->{$id} = \@results;
    }
    # Return the result.
    return $retVal;
}


=head3 contig_sequences

    my $contigHash = $sapObject->contig_sequences({
                        -ids => [$contig1, $contig2, ...]
                    });

Return the DNA sequences for the specified contigs.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of contig IDs. Note that the contig ID contains the
genome ID as a prefix (e.g. C<100226.1:NC_003888>).

=back

=item RETURN

Returns a reference to a hash that maps each contig ID to its DNA sequence.

    $contigHash = { $contig1 => $dna1, $contig2 => $dna2, ... };

=back

=cut

sub contig_sequences {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of contig IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs, getting the DNA sequences and putting them in the
    # return hash.
    for my $id (@$ids) {
        my @dna = $sap->GetFlat("Contig HasSection DNASequence",
                                'Contig(id) = ? ORDER BY DNASequence(id)', [$id],
                                'DNASequence(sequence)');
        $retVal->{$id} = join("", @dna);
    }
    # Return the result.
    return $retVal;
}


=head3 contig_lengths

    my $contigHash = $sapObject->contig_lengths({
                        -ids => [$contig1, $contig2, ...]
                    });

Return the lengths for the specified contigs.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of contig IDs. Note that the contig ID contains the
genome ID as a prefix (e.g. C<100226.1:NC_003888>).

=back

=item RETURN

Returns a reference to a hash that maps each contig ID to its length in base
pairs.

    $contigHash = { $contig1 => $len1, $contig2 => $len2, ... };

=back

=cut

sub contig_lengths {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of contig IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs, getting conting lengths and putting them in the
    # return hash.
    for my $id (@$ids) {
        my ($len) = $sap->GetFlat("Contig", 'Contig(id) = ?', [$id], 'length');
        $retVal->{$id} = $len;
    }
    # Return the result.
    return $retVal;
}

=head3 gene_correspondence_map

    my $geneHash =          $sapObject->gene_correspondence_map({
                                -genome1 => $genome1,
                                -genome2 => $genome2,
                                -fullOutput => 1,
                                -passive => 0
                            });

Return a map of genes in the specified second genome that correspond to genes in
the specified first genome.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -genome1

ID of the first genome of interest.

=item -genome2

ID of the second genome of interest.

=item -fullOutput (optional)

If C<1>, then instead of a simple hash map, a list of lists will be returned.
If C<2>, then the list will contain unidirectional correspondences from the target
back to the source as well as bidirectional corresopndences and unidirectional
correspondences from the source to the target. The default is C<0>, which returns
the hash map.

=item -passive (optional)

If TRUE, then an undefined value will be returned if no correspondence file
exists. If FALSE, a correspondence file will be created and cached on the server
if one does not already exist. This is an expensive operation, so set the flag
to TRUE if you are worried about performance. The default is FALSE.

=back

=item RETURN

This method will return an undefined value if either of the genome IDs is missing,
not found, or incomplete.

=over 8

=item Normal Output

Returns a hash that maps each gene in the first genome to a corresponding gene in
the second genome. The correspondence is determined by examining factors such as
functional role, conserved neighborhood, and similarity.

    $geneHash = { $g1gene1 => $g2gene1, $g1gene2 => $g2gene2,
                  $g1gene3 => $g2gene3, ... };

=item Output with -fullOutput >= 1

Returns a reference to list of sub-lists. Each sub-list contains 18 data items, as
detailed in L<ServerThing/Gene Correspondence List>.

=back

=back

=cut

sub gene_correspondence_map {
    # Get the parameters.
    my ($self, $args) = @_;
    # We'll put the results in here.
    my $retVal;
    # Get the two genome IDs.
    my $genome1 = $args->{-genome1};
    my $genome2 = $args->{-genome2};
    if (! defined $genome1) {
        Trace("-genome1 missing in gene_correspondence_map call.") if T(Corr => 1);
    } elsif (! defined $genome2) {
        Trace("-genome2 missing in gene_correspondence_map call.") if T(Corr => 1);
    } else {
        # We have genome IDs. Get the sapling database.
        my $sap = $self->{db};
        # Validate the genome IDs.
        my %completeMap =
            map { $_->[0] => $_->[1] } $sap->GetAll("Genome", "Genome(id) IN (?,?)",
                                                    [$genome1, $genome2], "id complete");
        if (! $completeMap{$genome1}) {
            Trace("Genome $genome1 not found or incomplete.") if T(Corr => 1);
        } elsif (! $completeMap{$genome2}) {
            Trace("Genome $genome2 not found or incomplete.") if T(Corr => 1);
        } else {
            # The genomes are both complete. Determine the output mode.
            my $fullOutput = $args->{-fullOutput} || 0;
            # Determine whether or not we're passive.
            my $passive = $args->{-passive};
            # This will hold the correspondence data.
            my $corrList = ServerThing::GetCorrespondenceData($genome1, $genome2, $passive,
                                                              $fullOutput == 2);
            # Do we have a result?
            if (defined $corrList) {
                # Check the output mode.
                if ($fullOutput) {
                    # Full output is the correspondence list itself.
                    $retVal = $corrList;
                } else {
                    # Normal output is a hash.
                    my %corrHash = map { $_->[0] => $_->[1] } @$corrList;
                    $retVal = \%corrHash;
                }
            } elsif (! $passive) {
                # Here we couldn't find a file AND the user is NOT in passive mode. That
                # indicates an error condition.
                Confess("Could not generate corresopndences from $genome1 to $genome2.");
            }
        }
    }
    # Return the result.
    return $retVal;
}


=head3 genome_contig_md5s

    my $genomeHash =        $sapObject->genome_contig_md5s({
                                -ids => [$genome1, $genome2, ...]
                            });

For each incoming genome, return a hash mapping its contigs to their MD5 identifiers.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the genome IDs.

=back

=item RETURN

Returns a hash that maps each incoming genome ID to a sub-hash that maps its contig IDs
to their MD5 identifiers. The MD5 identifiers are computed directly from the contig
DNA sequences.

    $genomeHash = { $genome1 => {$contig1a => $md5id1a, $contig1b => $md5id1b, ... },
                    $genome2 => {$contig2a => $md5id2a, $contig2b => $md5id2b, ... },
                    ... };

=back

=cut

sub genome_contig_md5s {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs, asking for contigs.
    for my $id (@$ids) {
        my @contigs = $sap->GetAll("IsMadeUpOf Contig", "IsMadeUpOf(from-link) = ?", [$id],
                                    ['to-link', 'Contig(md5-identifier)']);
        # If we found contigs for this genome, store them in the return hash. We must
        # convert the list of 2-tuples to a hash.
        if (@contigs) {
            $retVal->{$id} = { map { $_->[0] => $_->[1] } @contigs };
        }
    }
    # Return the result.
    return $retVal;
}

=head3 genome_contigs

    my $genomeHash =        $sapObject->genome_contigs({
                                -ids => [$genome1, $genome2, ...]
                            });

For each incoming genome, return a list of its contigs.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the genome IDs.

=back

=item RETURN

Returns a hash that maps each incoming genome ID to a list of its contig IDs.

    $genomeHash = { $genome1 => [$contig1a, $contig1b, ...],
                    $genome2 => [$contig2a, $contig2b, ...],
                    ... };

=back

=cut

sub genome_contigs {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs, asking for contigs.
    for my $id (@$ids) {
        my @contigs = $sap->GetFlat("IsMadeUpOf", "IsMadeUpOf(from-link) = ?", [$id],
                                    'to-link');
        # If we found contigs for this genome, store them in the return hash.
        if (@contigs) {
            $retVal->{$id} = \@contigs;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 genome_data

    my $genomeHash =        $sapObject->genome_data({
                                -ids => [$genome1, $genome2, ...],
                                -data => [$fieldA, $fieldB, ...]
                            });

Return the specified data items for the specified genomes.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of genome IDs.

=item -data

Reference to a list of data field names. The possible data field names are
given below.

=over 12

=item complete

C<1> if the genome is more or less complete, else C<0>.

=item contigs

The number of contigs for the genome

=item dna-size

The number of base pairs in the genome

=item domain

The domain of the genome (Archaea, Bacteria, ...).

=item gc-content

The amount of GC base pairs in the genome, expressed as a percentage of the
genome's DNA.

=item genetic-code

The genetic code used by this genome.

=item pegs

The number of protein encoding genes in the genome.

=item rnas

The number of RNAs in the genome.

=item name

The scientific name of the genome.

=item taxonomy

The genome's full taxonomy as a comma-separated string.

=item md5

The MD5 identifier computed from the genome's DNA sequences.

=back

=back

=item RETURN

Returns a hash mapping each incoming genome ID to an n-tuple. Each tuple
will contain the specified data fields for the computed gene in the specified
order.

    $genomeHash =  { $id1 => [$data1A, $data1B, ...],
                     $id2 => [$data2A, $data2B, ...],
                     ... };

=back

=cut

use constant GENOME_FIELDS => { complete => 'complete',
                                contigs => 'contigs',
                                'dna-size' => 'dna-size',
                                domain => 'domain',
                                'gc-content' => 'gc-content',
                                'genetic-code' => 'genetic-code',
                                pegs => 'pegs',
                                rnas => 'rnas',
                                name => 'scientific-name',
                                md5 => 'md5-identifier'};

use constant SPECIAL_FIELDS => { taxonomy => 1 };

sub genome_data {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of fields and perform a basic validation so we know we have a
    # list of data items.
    my $fields = $args->{-data};
    Confess("No data fields specified in \"genome_data\".") if ! defined $fields;
    Confess("Invalid data field list in \"genome_data\".") if ref $fields ne 'ARRAY';
    # There are two types of fields: GENOME_FIELDS contains the ones that
    # are actual fields in the genome record; SPECIAL_FIELDS are fields that
    # require special queries. The following hash maps each normal field's
    # database name to its output position.
    my %fieldNames;
    # This one maps each special field to its output position.
    my %otherFields;
    # Analyze the data field list and populate the field name hashes.
    for (my $i = 0; $i < @$fields; $i++) {
        # Get the current field's external name.
        my $field = $fields->[$i];
        # Check to see if this is a normal field.
        my $fieldName = GENOME_FIELDS->{$field};
        if (defined $fieldName) {
            # Yes. Store its database name in the field name map.
            $fieldNames{$fieldName} = $i;
        } else {
            # Check to see if this is a special field.
            my $found = SPECIAL_FIELDS->{$field};
            if ($found) {
                # Yes. Store it in the special name map.
                $otherFields{$field} = $i;
            } else {
                # No. It's a bad field.
                Confess("Invalid data field name \"$field\" in \"genome_data\".");
            }
        }
    }
    # Compute the list of normal field names.
    my @fieldNames = keys %fieldNames;
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the genomes.
    for my $id (@$ids) {
        # Get the normal data for this genome.
        my @tuple = $sap->GetEntityValues(Genome => $id, \@fieldNames);
        # Fill in the results.
        my @result;
        for (my $i = 0; $i < @fieldNames; $i++) {
            my $fieldName = $fieldNames[$i];
            $result[$fieldNames{$fieldName}] = $tuple[$i];
        }
        # Now run through and process the special fields.
        for my $otherField (keys %otherFields) {
            # Compute this field's value.
            if ($otherField eq 'taxonomy') {
                my @taxonomy = $sap->Taxonomy($id, 'names');
                $result[$otherFields{$otherField}] = join(",", @taxonomy);
            }
        }
        # Store the results.
        $retVal->{$id} = \@result;
    }
    # Return the result.
    return $retVal;
}


=head3 genome_domain

    my $genomeHash =        $sapObject->genome_domain({
                                -ids => [$genome1, $genome2, ...]
                            });

Return the domain for each specified genome (e.g. C<Archaea>, C<Bacteria>, C<Plasmid>).

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the genome IDs.

=back

=item RETURN

Returnss a hash that maps each incoming genome ID to its taxonomic domain.

=back

=cut

sub genome_domain {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return hash.
    my $retVal = {};
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs.
    for my $id (@$ids) {
        # Get the domain for this genome.
        my ($domain) = $sap->GetFlat('Genome', 'Genome(id) = ?', [$id], 'domain');
        # If we found it, store it in the return hash.
        if (defined $domain) {
            $retVal->{$id} = $domain;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 genome_fid_md5s

    my $genomeHash =        $sapObject->genome_fid_md5s({
                                -ids => [$genome1, $genome2, ...]
                            });

For each incoming genome, return a hash mapping its genes to their MD5 identifiers.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the genome IDs.

=back

=item RETURN

Returns a hash that maps each incoming genome ID to a sub-hash that maps its FIG feature IDs
to their MD5 identifiers. The MD5 identifiers are computed from the genome's MD5 identifier
and the gene's location in the genome.

    $genomeHash = { $genome1 => {$fid1a => $md5id1a, $fid1b => $md5id1b, ... },
                    $genome2 => {$fid2a => $md5id2a, $fid2b => $md5id2b, ... },
                    ... };

=back

=cut

sub genome_fid_md5s {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs, asking for contigs.
    for my $id (@$ids) {
        my @fids = $sap->GetAll("IsIdentifiedBy Identifier",
                                "IsIdentifiedBy(from-link) LIKE ? AND IsIdentifiedBy(to-link) LIKE ?",
                                ["fig|$id.%", "md5g|%"],
                                ['from-link', 'Identifier(natural-form)']);
        # If we found contigs for this genome, store them in the return hash. We must
        # convert the list of 2-tuples to a hash.
        if (@fids) {
            $retVal->{$id} = { map { $_->[0] => $_->[1] } @fids };
        }
    }
    # Return the result.
    return $retVal;
}


=head3 genome_ids

    my $genomeHash =        $sapObject->genome_ids({
                                -names => [$name1, $name2, ...],
                                -taxons => [$tax1, $tax2, ...]
                            });

Find the specific genome ID for each specified genome name or taxonomic number.
This method helps to find the correct version of a given genome when only the
species and strain are known.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -names (optional)

Reference to a list of genome scientific names, including genus, species, and
strain (e.g. C<Streptomyces coelicolor A3(2)>). A genome ID will be found (if any)
for each specified name.

=item taxons (optional)

Reference to a list of genome taxonomic numbers. These are essentially genome IDs
without an associated version number (e.g. C<100226>). A specific matching
genome ID will be found; the one chosen will be the one with the highest version
number that is not a plasmid.

=back

=item RETURN

Returns a hash mapping each incoming name or taxonomic number to the corresponding
genome ID.

    $genomeHash = { $name1 => $genome1, $name2 => $genome2, ...
                    $tax1 => $genome3, $tax2 => $genome4, ... };

=back

=cut

sub genome_ids {
    # Get the parameters.
    my ($self, $args) = @_;
    # We'll put our results in here.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of names.
    my $names = ServerThing::GetIdList(-names => $args, 1);
    # Loop through the names.
    for my $name (@$names) {
        # Find a genome ID for this name.
        my ($genome) = $sap->GetFlat('Genome', 'Genome(scientific-name) = ?', [$name],
                                     'id');
        # If we found one, store it in the return hash.
        if ($genome) {
            $retVal->{$name} = $genome;
        }
    }
    # Get the list of taxonomic IDs.
    my $taxons = ServerThing::GetIdList(-taxons => $args, 1);
    # Loop through the taxons.
    for my $taxon (@$taxons) {
        # Find the genome IDs for this taxon. Note we exclude plasmids.
        my (@genomes) = $sap->GetFlat('Genome',
                                      'Genome(id) LIKE ? AND Genome(domain) <> ?',
                                      ["$taxon.%", 'Plasmid'], 'id');
        # Only proceed if we found something.
        if (@genomes) {
            # Find the highest version number.
            my @sorted = sort { $a <=> $b } map { $_ =~ /\d+\.(.+)/; $1 } @genomes;
            # Use it to form the result.
            $retVal->{$taxon} = "$taxon.$sorted[$#sorted]";
        }
    }
    # Return all the results found.
    return $retVal;
}


=head3 genome_metrics

    my $genomeHash =        $sapObject->genome_metrics({
                                -ids => [$genome1, $genome2, ...]
                            });

For each incoming genome ID, returns the number of contigs, the total
number of base pairs in the genome's DNA, and the genome's default genetic
code.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of genome IDs.

=back

=item RETURN

Returns a hash mapping each incoming genome ID to a 3-tuple consisting of (0)
the number of contigs, (1) the total DNA size, and (2) the genome's default
genetic code.

    $genomeHash = { $genome1 => [$contigCount1, $baseCount1, $geneticCode1],
                    $genome2 => [$contigCount2, $baseCount2, $geneticCode2],
                    ... };

=back

=cut

sub genome_metrics {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs, getting the desired values.
    for my $id (@$ids) {
        my ($contigs, $dnaSize, $code) = $sap->GetEntityValues(Genome =>  $id,
                                                               [qw(contigs
                                                                dna-size
                                                                genetic-code)]);
        # Only proceed if we found this genome.
        if (defined $contigs) {
            $retVal->{$id} = [$contigs, $dnaSize, $code];
        }
    }
    # Return the result.
    return $retVal;
}

=head3 genome_names

    my $idHash =            $sapObject->genome_names({
                                -ids => [$id1, $id2, ...],
                                -numbers => 1
                            });

Return the name of the genome containing each specified feature or genome.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of identifiers. Each identifier can be a prefixed feature ID
(e.g. C<fig|100226.1.peg.3361>, C<uni|P0AC98>) or a genome ID (C<83333.1>,
C<360108.3>).

=item -numbers (optional)

If TRUE, the genome ID number will be returned instead of the name. Note that
this facility is only useful when the incoming identifiers are feature IDs,
as genome IDs would be mapped to themselves.

=back

=item RETURN

Returns a reference to a hash mapping each incoming feature ID to the scientific
name of its parent genome. If an ID refers to more than one real feature, only
the first feature's genome is returned.

    $idHash = { $id1 => $genomeName1, $id2 => $genomeName2, ... };

=back

=cut

sub genome_names {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Compute the output type (genome ID vs. name).
    my $outField = ($args->{-numbers} ? 'id' : 'scientific-name');
    # Get the sapling database object.
    my $sap = $self->{db};
    # Get the list of identifiers.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through them.
    for my $id (@$ids) {
        # The name we find will be put in here.
        my $name;
        # Is this a genome ID?
        if ($id =~ /^\d+\.\d+$/) {
            # Yes. Query the desired field from the genome record. Note that
            # if we are in number mode, this process will return the incoming
            # ID if the genome exists and an undefined value otherwise, which
            # is a meaningful result.
            ($name) = $sap->GetFlat("Genome", "Genome(id) = ?", [$id],
                                    "Genome($outField)");
        } else {
            # This is a feature identifier.
            ($name) = $sap->GetFlat("Identifies IsOwnedBy Genome",
                                    'Identifies(from-link) = ?', [$id],
                                    "Genome($outField)");
        }
        # If we found a result, store it in the return hash.
        if (defined $name) {
            $retVal->{$id} = $name;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 genomes_by_md5

    my $md5Hash =           $sapObject->genomes_by_md5({
                                -ids => [$md5id1, $md5id2, ...],
                                -names => 1
                            });

Find the genomes associated with each specified MD5 genome identifier. The MD5
genome identifier is computed from the DNA sequences of the genome's contigs; as
a result, two genomes with identical sequences arranged in identical contigs
will have the same MD5 identifier even if they have different genome IDs.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of MD5 genome identifiers.

=item -names (optional)

If TRUE, then both genome IDs and their associated names will be returned;
otherwise, only the genome IDs will be returned. The default is FALSE.

=back

=item RETURN

Returns a reference to a hash keyed by incoming MD5 identifier. Each identifier
maps to a list of genomes. If C<-names> is FALSE, then the list is of genome IDs;
if C<-names> is TRUE, then the list is of 2-tuples, each consisting of (0) a genome
ID and (1) the associated genome's scientific name.

=over 8

=item if C<-names> = TRUE

    $md5Hash = { $md5id1 => [[$genome1a, $name1a], [$genome1b, $name1b], ...],
                 $md5id2 => [[$genome2a, $name2a], [$genome2b, $name2b], ...],
                 ... };

=item if C<-names> = FALSE

    $md5Hash = { $md5id1 => [$genome1a, $genome1b, ...],
                 $md5id2 => [$genome2a, $genome2b, ...],
                 ... };

=back

=back

=cut

sub genomes_by_md5 {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return hash.
    my $retVal = {};
    # Get the list of incoming IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Determine whether or not we want the scientific names.
    my $namesFlag = $args->{-names} || 0;
    # Compute the field list depending on whether or not we want the scientific
    # names.
    my @fields = 'id';
    if ($namesFlag) {
        push @fields, 'scientific-name';
    }
    # Loop through the incoming MD5 ids.
    for my $id (@$ids) {
        # Get the genomes for this ID.
        my @rows = $sap->GetAll("Genome", 'Genome(md5-identifier) = ?', [$id], \@fields);
        # Store the results depending on the mode.
        if ($namesFlag) {
            # When we're asking for names, this is easy, because the output is in
            # exactly the form we want.
            $retVal->{$id} = \@rows
        } else {
            # When we only want IDs, we have to dereference the sub-lists so that
            # we have a list of strings instead of a list of 1-tuples.
            $retVal->{$id} = [map { $_->[0] } @rows];
        }
    }
    # Return the results.
    return $retVal;
}

=head3 intergenic_regions

    my $locList =           $sapObject->intergenic_regions({
                                -genome => $genome1,
                                -type => ['peg', 'rna']
                            });

Return a list of L</Location Strings> for the regions in the specified genome that are
not occupied by genes of the specified types. All of these will be construed to be on
the forward strand, and sorted by contig ID and start location within contig.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -genome

ID of the genome whose intergenic regions are to be returned.

=item -type (optional)

Reference to a list of gene types. Only genes of the specified type will be considered
to be occupying space on the contigs. Typically, this parameter will either be C<peg>
or a list consisting of C<peg> and C<rna>. The default is to allow all gene types,
but this will not generally produce a good result.

=back

=item RETURN

Returns a reference to a list of location strings, indicating the intergenic region
locations for the genome.

    $locList = [$loc1, $loc2, ...]

=back

=cut

sub intergenic_regions {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # The regions found will be stored in here.
    my @retVal;
    # Get the genome ID.
    my $genome = $args->{-genome};
    Confess("No genome specified for intergenic_regions.") if ! defined $genome;
    # Get the list of gene types. If none are specified, we'll get an empty list.
    my $types = ServerThing::GetIdList(-type => $args, 1);
    # We need to create the filter for the feature lookup. This list will contain
    # all of the parameters except the contig ID.
    my @filterParms;
    # Each parameter will be a feature ID pattern that only matches genes of the
    # desired type for our genome.
    for my $type (@$types) {
        push @filterParms, "fig|$genome.$type.%";
    }
    # Start the filter string with the contig ID.
    my $filterString = "IsLocusFor(from-link) = ?";
    if (scalar @filterParms) {
        # Here we have additional filtering by feature type.
        my @clauses;
        for my $filterParm (@filterParms) {
            push @clauses, "IsLocusFor(to-link) LIKE ?"
        }
        $filterString .= " AND (" . join(" OR ", @clauses) . ")";
    }
    # Finish off with an ordering.
    $filterString .= " ORDER BY IsLocusFor(from-link), IsLocusFor(begin)";
    # Now we have everything we need to create queries for the occupied regions of
    # a contig. The next step is to get the contigs.
    my @contigs = $sap->GetAll("IsMadeUpOf Contig", "IsMadeUpOf(from-link) = ?",
                                [$genome], 'Contig(id) Contig(length)');
    # Loop through the contigs.
    for my $contigData (@contigs) {
        # Get the contig ID and length.
        my ($contig, $contigLen) = @$contigData;
        # Denote that our current position on the contig is the first base pair.
        my $loc = 1;
        # Create a query to get all the occupied regions of the contig.
        my $query = $sap->Get("IsLocusFor", $filterString, [$contig, @filterParms]);
        # Loop through the results.
        while (my $region = $query->Fetch()) {
            # Get the start and length of this occupied region.
            my $begin = $region->PrimaryValue('begin');
            my $len = $region->PrimaryValue('len');
            # Is there an intergenic region before the start of this new area?
            if ($begin > $loc) {
                # Yes, write it out.
                my $regionLen = $begin - $loc;
                push @retVal, $contig . "_$loc+$regionLen";
            }
            # Record the end of this region as the last occupied position, if it's past
            # the current limit.
            my $regionLast = $begin + $len;
            $loc = $regionLast if $loc < $regionLast;
        }
        # Check for a residual at the end of the contig.
        if ($contigLen > $loc) {
            my $regionLen = $contigLen - $loc;
            push @retVal, $contig . "_$loc+$regionLen";
        }
    }
    # Return the list of locations found.
    return \@retVal;
}

=head3 is_prokaryotic

    my $genomeHash =        $sapObject->is_prokaryotic({
                                -ids => [$genome1, $genome2, ...]
                            });

For each incoming genome ID, returns 1 if it is prokaryotic and 0
otherwise.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the relevant genome IDs.

=back

=item RETURN

Returns a reference to a hash that maps each incoming genome ID to C<1> if it is
a prokaryotic genome and C<0> otherwise.

    $genomeHash = { $genome1 => $flag1, $genome2 => $flag2, ... };

=back

=cut

sub is_prokaryotic {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the incoming IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the IDs, pulling the prokaryotic flag.
    for my $id (@$ids) {
        my ($flag) = $sap->GetFlat("Genome", "Genome(id) = ?", [$id], 'Genome(prokaryotic)');
        $retVal->{$id} = $flag;
    }
    # Return the result.
    return $retVal;
}


=head3 mapped_genomes

    my $genomeHash =        $sapObject->mapped_genomes({
                                -ids => [$genome1, $genome2, ...]
                            });

For each incoming genome, return a list of the genomes that have an existing
gene correspondence map (see L<ServerThing/Gene Correspondence List>). Gene
correspondence maps indicate which genes in the target genome are the best hit
of each gene in the source genome. If a correspondence map does not yet exist,
it will be created when you ask for it, but this is an expensive process and it
is sometimes useful to find an alternate genome that will give you a faster
result.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the IDs for the genomes of interest. A (possibly
empty) result list will be returned for each one.

=back

=item RETURN

Returns a reference to a hash mapping each incoming genome ID to a list of
the IDs for the genomes which have existing correspondence maps on the
server.

    $genomeHash = { $genome1 => [$genome1a, $genome1b, ...],
                    $genome2 => [$genome2a, $genome2b, ...],
                    ... };

=back

=cut

sub mapped_genomes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get a hash of complete genomes. This is used to filter old, obsolete
    # genome IDs out of the output and to determine which directories we want
    # to examine.
    my $genomeHash = $self->all_genomes({ -complete => 1 });
    # Get the list of incoming genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the list.
    for my $id (@$ids) {
        # Get all of the genomes in this genome's correspondence map directory.
        my $orgDir = ServerThing::ComputeCorrespondenceDirectory($id);
        my @files = grep { exists $genomeHash->{$_} } OpenDir($orgDir, 0, 1);
        # Put them in the result list.
        $retVal->{$id} = \@files;
        # The correspondence maps are reversible, so we only keep half of them.
        # Our next task, then, is to find the converse correspondence files in
        # the directories of other genomes. The other genomes will be in directories
        # for genomes that satisfy the "must-flip" criterion when placed next to
        # this one.
        for my $otherID (keys %$genomeHash) {
            if (ServerThing::MustFlipGenomeIDs($id, $otherID)) {
                # Here we have an ID for a corresponding genome that will have us
                # in its directory instead of being in our directory.
                my ($fileName) = ServerThing::ComputeCorrespondenceFileName($id, $otherID);
                if (-f $fileName) {
                    push @{$retVal->{$id}}, $otherID;
                }
            }
        }
    }
    # Return the result hash.
    return $retVal;
}

=head3 otu_members

    my $genomeHash =        $sapObject->otu_members({
                                -ids => [$genome1, $genome2, ...]
                            });

For each incoming genome, return the name and ID of each other genome in the same
OTU.

=over 4

=item parameter

The parameter shoudl be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the IDs for the genomes of interest.

=back

=item RETURN

Returns a reference to a hash mapping each incoming genome ID to a sub-hash.
The sub-hash is keyed by genome ID, and maps the ID of each genome in the same
OTU to its name.

    $genomeHash = { $genome1 => { $genome1a => $name1a, $genome1b => $name1b, ... },
                    $genome2 => { $genome2a => $name2a, $genome2b => $name2b, ... },
                    ... };

=back

=cut

sub otu_members {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return hash.
    my $retVal = {};
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the genomes.
    for my $genomeID (@$ids) {
        # Get the related genomes.
        my %neighbors = map { $_->[0] => $_->[1] }
                    $sap->GetAll("IsCollectedInto IsCollectionOf Genome",
                        'IsCollectedInto(from-link) = ? AND Genome(id) <> ?',
                        [$genomeID,$genomeID], [qw(Genome(id) Genome(scientific-name))]);
        # Store them in the return hash.
        $retVal->{$genomeID} = \%neighbors;
    }
    # Return the result.
    return $retVal;
}

=head3 representative

    my $genomeHash =        $sapObject->representative({
                                -ids => [$genome1, $genome2, ...]
                            });

Return the representative genome for each specified incoming genome ID.
Genomes with the same representative are considered closely related, while
genomes with a different representative would be considered different
enough that similarities between them have evolutionary significance.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the IDs for the genomes of interest.

=back

=item RETURN

Returns a reference to a hash mapping each incoming genome ID to the ID of
its representative genome.

    $genomeHash = { $genome1 => $genome1R, $genome2 => $genome2R, ... };

=back

=cut

sub representative {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database object.
    my $sap = $self->{db};
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the genomes.
    for my $genomeID (@$ids) {
        # Get this genome's representative.
        my ($representative) = $sap->GetFlat("IsCollectedInto IsCollectionOf",
                                            'IsCollectedInto(from-link) = ? AND IsCollectionOf(representative) = ?',
                                            [$genomeID,'1'],
                                            'IsCollectionOf(to-link)');
        # Only proceed if we found one. If we didn't, it means the genome ID
        # is invalid.
        if ($representative) {
            $retVal->{$genomeID} = $representative;
        }
    }
    # Return the result.
    return $retVal;
}


=head3 representative_genomes

    my $mappings =          $sapObject->representative_genomes();

Compute mappings for the genome sets (OTUs) in the database. This method will
return a mapping from each genome to its genome set ID and from each
genome set ID to a list of the genomes in the set. For the second
mapping, the first genome in the set will be the representative.

This method does not require any parameters.

=over 4

=item RETURN

Returns a reference to a 2-tuple. The first element is a reference to a hash
mapping genome IDs to genome set IDs; the second element is a reference to a hash
mapping each genome set ID to a list of the genomes in the set. The first genome
in each of these lists will be the set's representative.

    $mappings = [ { $genome1 => $set1, $genome2 => $set2, ... },
                  { $set1 => [$genome1R, $genome1a, $genome1b, ...],
                    $set2 => [$genome2R, $genome2a, $genome2b, ...],
                    ... }
                ];

=back

=cut

sub representative_genomes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variables.
    my %genomes_to_otus;
    my %otus_to_genomes;
    # Get the sapling database.
    my $sap = $self->{db};
    # Read in the genome sets. The ordering insures that for each set, we see the
    # representative genome first. This means we don't need to do any fancy
    # checking when we build the list of genomes in each set.
    my @setData = $sap->GetAll("IsCollectionOf Genome",
                               'ORDER BY IsCollectionOf(from-link), IsCollectionOf(representative) DESC, IsCollectionOf(to-link)',
                               [], [qw(from-link to-link)]);
    # Loop through the set data returned.
    for my $setDatum (@setData) {
        # Get the genome data.
        my ($setID, $genome) = @$setDatum;
        # Add this genome to the two hashes.
        push @{$otus_to_genomes{$setID}}, $genome;
        $genomes_to_otus{$genome} = $setID;
    }
    # Return the result.
    return [\%genomes_to_otus, \%otus_to_genomes];
}

=head3 submit_gene_correspondence

    my $statusCode =    $sapObject->submit_gene_correspondence({
                            -genome1 => $genome1,
                            -genome2 => $genome2,
                            -correspondences => $corrList,
                            -passive => 1
                        });

Submit a set of gene correspondences to be stored on the server.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -genome1

ID of the source genome for the correspondence.

=item -genome2

ID of the target genome for the correspondence.

=item -correspondences

Reference to a list of lists containing the correspondence data
(see L<ServerThing/Gene Correspondence List>).

=item -passive (optional)

If TRUE, then the file will not be stored if one already exists. If FALSE, an
existing correspondence file will be overwritten. The default is FALSE.

=back

=item RETURN

Returns TRUE (C<1>) if the correspondences were successfully stored, FALSE
(C<0>) if they were rejected or an error occurred.

=back

=cut

sub submit_gene_correspondence {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable. We assume an error unless successful.
    my $retVal = 0;
    # Get the two genome IDs.
    my $genome1 = $args->{-genome1};
    if (! defined $genome1) {
        Confess("-genome1 missing in submit_gene_correspondence call.");
    }
    my $genome2 = $args->{-genome2};
    if (! defined $genome2) {
        Confess("-genome2 missing in submit_gene_correspondence call.");
    }
    # Get the correspondence list.
    my $corrList = $args->{-correspondences};
    if (! defined $corrList) {
        Confess("-correspondences missing in submit_gene_correspondence call.");
    }
    # Compute the name for the genome correspondence file and find out if the
    # genome IDs are in the right order.
    my ($fileName, $genomeA, $genomeB) = ServerThing::ComputeCorrespondenceFileName($genome1, $genome2);
    my $converse = ($genomeA ne $genome1);
    # Determine if we're active or passive. If we're passive and the file already
    # exists, we simply return success and quit.
    my ($existingFileName) = ServerThing::CheckForGeneCorrespondenceFile($genome1, $genome2);
    if ($args->{-passive} && $existingFileName) {
        $retVal = 1;
        Trace("Correspondence for $genome1 to $genome2 already exists. Skipped in passive mode.") if T(Corr => 3);
    } else {
        # Insure the correspondence list is valid.
        if (ref $corrList ne 'ARRAY') {
            Trace("Invalid correspondence list in submit_gene_correspondence for $genome1 to $genome2: not an array.") if T(Corr => 0);
        } else {
            # Loop through the list, checking for errors.
            my $errorCount = 0;
            for (my $i = 0; $i < scalar(@$corrList); $i++) { last if $errorCount > 0;
                my $row = $corrList->[$i];
                if (ref $row ne 'ARRAY') {
                    Trace("Invalid correspondence list in submit_gene_correspondence for $genome1 to $genome2: row $i is not an array.") if T(Corr => 0);
                    $errorCount++;
                } else {
                    $errorCount += ServerThing::ValidateGeneCorrespondenceRow($row);
                    if ($errorCount) {
                        Trace("Invalid correspondence list in submit_gene_correspondence for $genome1 to $genome2: row $i has errors.") if T(Corr => 0);
                    } elsif ($converse) {
                        # Here we have to flip the row to get it in the right order.
                        ServerThing::ReverseGeneCorrespondenceRow($row);
                    }
                }
            }
            if (! $errorCount) {
                # Now we need to verify the genome IDs.
                my $sap = $self->{db};
                for my $genome ($genome1, $genome2) {
                    my ($complete) = $sap->GetFlat('Genome', 'Genome(id) = ?', [$genome],
                                                   'complete');
                    if (! $complete) {
                        Trace("$genome missing or incomplete. Cannot store correspondence file.") if T(Corr => 0);
                        $errorCount++;
                    }
                }
                if (! $errorCount) {
                    # Now we know we can store the correspondence data. Try to open a temporary
                    # file to hold the data,
                    my $tempFileName = "$fileName.$$.tmp";
                    my $oh;
                    if (! open($oh, ">$tempFileName")) {
                        Trace("Could not open correspondence temp file: $!") if T(Corr => 0);
                    } else {
                        # Store the data in the file.
                        for my $row (@$corrList) {
                            print $oh join("\t", @$row) . "\n";
                        }
                        # Close the temporary file.
                        if (! close $oh) {
                            Trace("Error closing $tempFileName. Correspondence store aborted.") if T(Corr => 0);
                        } else {
                            # Try to rename it.
                            if (rename $tempFileName, $fileName) {
                                # It worked! Fix the permissions and denote success.
                                chmod 0664, $fileName;
                                $retVal = 1;
                            } else {
                                Trace("Error renaming $tempFileName to $fileName. Correspondence store aborted.") if T(Corr => 0);
                            }
                        }
                    }
                    # Insure the temporary file is deleted.
                    if (-f $tempFileName) {
                        unlink $tempFileName;
                    }
                }
            }
        }
    }
    # Return the success indicator.
    return $retVal;
}


=head3 taxonomy_of

    my $genomeHash =    $sapObject->taxonomy_of({
                            -ids => [$genome1, $genome2, ...],
                            -format => 'numbers'
                        });

Return the taxonomy of each specified genome. The taxonomy will start at
the domain level and moving down to the node where the genome is
attached.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of genome IDs. A taxonomy will be generated for each
specified genome.

=item -format (optional)

Format for the elements of the taxonomy string. If C<numbers>, then each
taxonomy element will be represented by its number; if C<names>, then each
taxonomy element will be represented by its primary name; if C<both>, then
each taxonomy element will be represented by a number followed by the name.
The default is C<names>.

=back

=item RETURN

Returns a reference to a hash mapping incoming genome IDs to taxonomies.
Each taxonomy will be a list of strings, starting from the domain and
ending with the genome.

=over 8

=item Normal Output

    $genomeHash = { $genome1 => [$name1a, $name1b, ...],
                    $genome2 => [$name2a, $name2b, ...],
                    ... };

=item Output if -format = numbers

    $genomeHash = { $genome1 => [$num1a, $num1b, ...],
                    $genome2 => [$num2a, $num2b, ...],
                    ... };

=item Output if =format = both

    $genomeHash = { $genome1 => ["$num1a $name1a", "$num1b $name1b", ...],
                    $genome2 => ["$num2a $name2a", "$num2b $name2b", ...],
                    ... };

=back

=back

=cut

sub taxonomy_of {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database object.
    my $sap = $self->{db};
    # Determine the format.
    my $format = $args->{-format} || 'names';
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the genomes.
    for my $genomeID (@$ids) {
        # Get this genome's taxonomy and put it in the return hash.
        $retVal->{$genomeID} = [ $sap->Taxonomy($genomeID, $format) ];
    }
    # Return the result.
    return $retVal;
}


=head2 Scenario Data Methods

=head3 scenario_names

    my $scenarioHash =      $sapObject->scenario_names({
                                -subsystem => $subsys1
                            });

Return the names of all the scenarios for the specified subsystem. Each scenario
has an internal ID number and a common name. This method returns both.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -subsystem

Name of the subsystem whose scenarios are desired.

=back

=item RETURN

Returns a hash mapping the ID numbers of the subsystem's scenarios to their
common names.

    $scenarioHash = { $id1 => $name1, $id2 => $name2, ... };

=back

=cut

sub scenario_names {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my %retVal;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the ID of the target subsystem.
    my $subsystem = $args->{-subsystem};
    if (! defined $subsystem) {
        Confess("No -subsystem specified for scenario_names.");
    } else {
        # Retrieve the scenario names and put them into the result hash.
        %retVal = map { $_->[0] => $_->[1] }
                    $sap->GetAll("Subsystem IsSubInstanceOf Scenario",
                                 'Subsystem(id) = ?', [$subsystem], [qw(Scenario(id)
                                 Scenario(common-name))]);
    }
    # Return the result hash.
    return \%retVal;
}


=head2 Subsystem Data Methods

=head3 all_subsystems

    my $subsysHash =        $sapObject->all_subsystems({
                                -usable => 1,
                                -exclude => [$type1, $type2, ...],
                                -aux => 1
                            });

Return a list of all subsystems in the system. For each subsystem, this
method will return the ID, curator, the classifications, and roles.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys, all of
which are optional. Because all of the keys are optional, it is permissible to
pass an empty hash or no parameters at all.

=over 8

=item -usable (optional)

If TRUE, then only results from usable subsystems will be included. If FALSE,
then results from all subsystems will be included. The default is TRUE.

=item -exclude (optional)

Reference to a list of special subsystem types that should be excluded from the
result list. The permissible types are C<cluster-based> and C<experimental>.
Normally cluster-based subsystems are included, but experimental subsystems
are only included if the C<-usable> option is turned off.

=item -aux (optional)

If TRUE, then auxiliary roles will be included in the output. The default is
FALSE, meaning they will be excluded.

=back

=item RETURN

Returns a hash mapping each subsystem ID to a 3-tuple consisting of (0) the name of the
curator, (1) a reference to a list of the subsystem classifications, and (2) a reference
to a list of the subsystem's roles.

    $subsysHash = { $sub1 => [$curator1, [$class1a, $class1b, ...], [$role1a, $role1b, ...]],
                    $sub2 => [$curator2, [$class2a, $class2b, ...], [$role2a, $role2b, ...]],
                    ... };

=back

=cut

sub all_subsystems {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the spaling database.
    my $sapling = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Compute the filter based on the parameters.
    my $filter = "";
    ServerThing::AddSubsystemFilter(\$filter, $args, 1);
    # Create a hash for walking up the subsystem class hierarchy.
    my %classMap = map { $_->[0] => $_->[1] } $sapling->GetAll("IsSubclassOf",
                                                               "", [],
                                                               [qw(from-link to-link)]);
    # Read the subsystem role data from the database.
    my @roleData = $sapling->GetAll("Subsystem Includes Role AND Subsystem IsInClass SubsystemClass",
                                    $filter, [],
                                    [qw(Subsystem(id) Subsystem(curator)
                                        SubsystemClass(id) Role(id))]);
    # Loop through the subsystems, building the result hash.
    for my $roleDatum (@roleData) {
        my ($subsystem, $curator, $class, $role) = @$roleDatum;
        # Is this subsystem new?
        if (! exists $retVal->{$subsystem}) {
            # Yes. Get its classification data. We trace the classifications from
            # the bottom up, so new ones are shifted onto the front.
            my @classes;
            while ($class) {
                unshift @classes, $class;
                $class = $classMap{$class};
            }
            # Create its hash entry.
            $retVal->{$subsystem} = [$curator, \@classes, []];
        }
        # Now we know an entry exists for this subsystem. Push this role onto it.
        push @{$retVal->{$subsystem}[2]}, $role;
    }
    # Return the result.
    return $retVal;
}

=head3 classification_of

    my $subsysHash =        $sapObject->classification_of({
                                -ids => [$sub1, $sub2, ...]
                            });

Return the classification for each specified subsystem.

=over 4

=item parameter

Reference to a hash of parameters with the following possible keys.

=over 8

=item -ids

Reference to a list of subsystem IDs.

=back

=item RETURN

Returns a hash mapping each incoming subsystem ID to a list reference. Each
list contains the classification names in order from the largest classification to
the most detailed.

    $subsysHash = { $sub1 => [$class1a, $class1b, ...],
                    $sub2 => [$class2a, $class2b, ...],
                    ... };

=back

=cut

sub classification_of {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of subsystem IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the subsystem IDs, getting the classification data.
    for my $id (@$ids) {
        # We'll build the classification list in here.
        my @classes;
        # Normalize the ID.
        my $realID = $sap->SubsystemID($id);
        # Get the low-level class.
        my ($class) = $sap->GetFlat("Subsystem IsInClass SubsystemClass",
                                    "Subsystem(id) = ?", [$realID], 'SubsystemClass(id)');
        # Loop through the remaining classes. Note that since we're moving up
        # the hierarchy, new classes are added at the beginning.
        while (defined $class) {
            unshift @classes, $class;
            ($class) = $sap->GetFlat("SubsystemClass IsSubclassOf SubsystemClass2",
                                     "SubsystemClass(id) = ?", [$class],
                                     'SubsystemClass2(id)');
        }
        # Store this classification.
        $retVal->{$id} = \@classes;
    }
    # Return the result.
    return $retVal;
}


=head3 genomes_to_subsystems

    my $genomeHash =        $sapObject->genomes_to_subsystems({
                                -ids => [$genome1, $genome2, ...],
                                -all => 1,
                                -usable => 0,
                                -exclude => ['cluster-based', 'experimental', ...]
                            });

Return a list of the subsystems participated in by each of the specified
genomes.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the genome IDs.

=item -all (optional)

If TRUE, all subsystems will be returned, including those in which the genome
does not appear to implement the subsystem and those in which the subsystem
implementation is incomplete. The default is FALSE, in which case only subsystems
that are completely implemented by the genome will be returned.

=item -usable (optional)

If TRUE, then only results from usable subsystems will be included. If FALSE,
then results from all subsystems will be included. The default is TRUE.

=item -exclude (optional)

Reference to a list of special subsystem types that should be excluded from the
result list. The permissible types are C<cluster-based> and C<experimental>.
Normally cluster-based subsystems are included, but experimental subsystems
are only included if the C<-usable> option is turned off.

=back

=item RETURN

Returns a hash mapping each genome ID to a list of 2-tuples. Each 2-tuple will
contain a subsystem name followed by a variant code.

    $genomeHash = { $genome1 => [[$sub1a, $variantCode1a], [$sub1b, $variantCode1b], ...],
                    $genome2 => [[$sub2a, $variantCode2a], [$sub2b, $variantCode2b], ...],
                    ... };

=back

=cut

sub genomes_to_subsystems {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of genome IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Get the all-subsystems flag.
    my $all = $args->{-all} || 0;
    # Format the filter clause. If we're getting all subsystems, it's shorter than
    # the norm.
    my $filterString = 'Genome(id) = ?';
    my @parms = ();
    if (! $all) {
        $filterString .= ' AND Variant(type) = ?';
        push @parms, 'normal';
    }
    # Add subsystem type filtering.
    ServerThing::AddSubsystemFilter(\$filterString, $args);
    # Loop through the genome IDs.
    for my $genomeID (@$ids) {
        # Get the subsystems for this genome.
        my @data = $sap->GetAll("Genome Uses Implements Variant IsDescribedBy Subsystem",
                                $filterString, [$genomeID, @parms],
                                [qw(Subsystem(id) Variant(code))]);
        # If we found any, put them in the result hash.
        if (scalar @data) {
            $retVal->{$genomeID} = \@data;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 get_subsystems

    my $subsysHash =        $sapObject->get_subsystems({
                                -ids => [$sub1, $sub2, ...]
                            });

Get a complete description of each specified subsystem. This will include
the basic subsystem properties, the list of roles, and the spreadsheet.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of subsystem IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming subsystem ID to a sub-hash
that completely describes the subsystem. The keys for the sub-hash are as
follows.

=over 8

=item curator

The name of the subsystem's curator.

=item version

The subsystem's current version number.

=item notes

The text of the subsystem notes.

=item desc

The description of the subsystem.

=item roles

Reference to a list of 3-tuples, one for each role in the subsystem. Each
3-tuple will contain (0) the role abbreviation, (1) C<1> if the role is
auxiliary and C<0> otherwise, and (2) the ID (name) of the role.

=item spreadsheet

Reference to a list of 5-tuples. For each molecular machine implementing the
subsystem, there is a 5-tuple containing (0) the target genome ID, (1) the
relevant region string, (2) C<1> if the molecular machine is curated and C<0>
if it was computer-assigned, (3) the variant code for the implemented variant,
and (4) a reference to a list of sub-lists, one per role (in order), with each
sub-list containing the IDs of all features performing that role.

=back

    $subsysHash = { $sub1 =>
                        { curator => $curator1,
                          version => $version1,
                          notes => $notes1,
                          desc => $desc1,
                          roles => [[$abbr1a, $aux1a, $role1a],
                                    [$abbr1b, $aux1b, $role1b], ... ],
                          spreadsheet => [
                            [$genome1x, $region1x, $curated1x, $variant1x,
                                [[$fid1xa1, $fid1xa2, ...], [$fid1xb1, $fid1xb2, ...], ...]],
                            [$genome1y, $region1y, $curated1y, $variant1y,
                                [[$fid1ya1, $fid1ya2, ...], [$fid1yb1, $fid1yb2, ...], ...]],
                            ... ]
                        },
                    $sub2 =>
                        { curator => $curator2,
                          version => $version2,
                          notes => $notes2,
                          desc => $desc2,
                          roles => [[$abbr2a, $aux2a, $role2a],
                                    [$abbr2b, $aux2b, $role2b], ... ],
                          spreadsheet => [
                            [$genome2x, $region2x, $curated2x, $variant2x,
                                [[$fid2xa1, $fid2xa2, ...], [$fid2xb1, $fid2xb2, ...], ...]],
                            [$genome2y, $region2y, $curated1y, $variant1y,
                                [[$fid1ya1, $fid1ya2, ...], [$fid1yb1, $fid1yb2, ...], ...]],
                            ... ]
                        },

=back

=cut

sub get_subsystems {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of subsystems.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Create the return hash.
    my $retVal = {};
    # Loop through the subsystems.
    for my $id (@$ids) {
        # Get the basic subsystem data.
        my ($subData) = $sap->GetAll("Subsystem", 'Subsystem(id) = ?', [$id],
                                     [qw(curator description notes version)]);
        # Only proceed if the subsystem exists.
        if ($subData) {
            # Create the subsystem's sub-hash.
            my %subHash = ( curator => $subData->[0], desc => $subData->[1],
                            notes => $subData->[2], version => $subData->[3] );
            # Get the role list.
            my @roleList = $sap->GetAll("Includes",
                        'Includes(from-link) = ? ORDER BY Includes(sequence)',
                        [$id], [qw(abbreviation auxiliary to-link)]);
            $subHash{roles} = \@roleList;
            # Now we need to get the spreadsheet. We start with a list of the
            # genomes and their variants.
            my @machines = $sap->GetAll("Describes Variant IsImplementedBy MolecularMachine IsUsedBy",
                                        'Describes(from-link) = ?', [$id],
                                        [qw(IsUsedBy(to-link) MolecularMachine(region)
                                            MolecularMachine(curated) Variant(code))]);
            # Get the MD5 of the subsystem ID.
            my $subsysMD5 = ERDB::DigestKey($id);
            # Get all of the features in the subsystem, and organize them into
            # a hash.
            my %cells;
            my $qh = $sap->Get("Contains", 'Contains(from-link) LIKE ?', [$subsysMD5 . '%']);
            while (my $resultRow = $qh->Fetch()) {
                my $fromLink = $resultRow->PrimaryValue('from-link');
                my $fid = $resultRow->PrimaryValue('to-link');
                my (undef, $genome, $region, $role) = split /:/, $fromLink;
                push @{$cells{"$genome:$region"}{$role}}, $fid;
            }
            # Get the list of role abbreviations.
            my @roles = map { $_->[0] } @roleList;
            # Loop through the machines. For each machine, we must create the
            # list of spreadsheet cells and add it to the end.
            for my $machine (@machines) {
                # Get the sub-hash for this machine.
                my $machineH = $cells{"$machine->[0]:$machine->[1]"};
                # Map the cells into a list.
                my @row = map { $machineH->{$_} || [] } @roles;
                # Add the list to the machine.
                push @$machine, \@row;
            }
            # Store the machines in the sub-hash.
            $subHash{spreadsheet} = \@machines;
            # Store the subhash in the result hash.
            $retVal->{$id} = \%subHash;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 ids_in_subsystems

    my $subsysHash =        $sapObject->ids_in_subsystems({
                                -subsystems => [$sub1, $sub2, ...],
                                -genome => $genome1,
                                -grouped => 1,
                                -roleForm => 1,
                                -source => 'UniProt'
                            });

Return the features of each specified subsystems in the specified genome, or
alternatively, return all features of each specified subsystem.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -subsystems

Reference to a list of the IDs for the desired subsystems.

=item -genome (optional)

ID of the relevant genome, or C<all> to return the genes in all genomes for
the subsystem. The default is C<all>.

=item -grouped (optional)

If specified, then instead of being represented in a list, the feature IDs will be
represented in a comma-delimited string.

=item -roleForm (optional)

If C<abbr>, then roles will be represented by the role abbreviation; if C<full>, then
the role will be represented by its full name; if C<none>, then roles will not be
included and there will only be a single level of hashing-- by subsystem ID. The
default is C<abbr>.

=item -source (optional)

Database source for the output IDs-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
The default is C<SEED>.

=back

=item RETURN

Returns a hash mapping each subsystem ID to a sub-hash. Each sub-hash
maps the roles of the subsystem to lists of feature IDs. The roles are recorded
as role abbreviations.

=over 8

=item Normal Output

    $subsysHash = { $sub1 => { $roleAbbr1A => [$fid1Ax, $fid1Ay, ...],
                               $roleAbbr1B => [$fid1Bx, $fid1By, ...],
                               ... },
                    $sub2 => { $roleAbbr2A => [$fid2Ax, $fid2Ay, ...],
                               $roleAbbr2B => [$fid2Bx, $fid2By, ...],
                               ... },
                    ... };

=item Output if -roleForm = full

    $subsysHash = { $sub1 => { $role1A => [$fid1Ax, $fid1Ay, ...],
                               $role1B => [$fid1Bx, $fid1By, ...],
                               ... },
                    $sub2 => { $role2A => [$fid2Ax, $fid2Ay, ...],
                               $role2B => [$fid2Bx, $fid2By, ...],
                               ... },
                    ... };

=item Output if -roleForm = none

    $subsysHash = { $sub1 => [$fid1a, $fid1b, ...],
                    $sub2 => [$fid2a, $fid2b, ...],
                    ... };

=back

=back

=cut

sub ids_in_subsystems {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the subsystem ID list.
    my $subsystems = ServerThing::GetIdList(-subsystems => $args);
    # Get the genome ID.
    my $genome = $args->{-genome} || 'all';
    # The following list forms the tail of the machine role ID pattern. If a
    # specific genome is specified, then the pattern ends with the genome ID and
    # a wild card; otherwise, it's just a wild card.
    my @patternTail;
    push @patternTail, $genome if ($genome ne 'all');
    push @patternTail, "%";
    # Process the flags and options.
    my $grouped = $args->{-grouped} || 0;
    my $roleForm = $args->{-roleForm} || 'abbr';
    my $source = $args->{-source} || 'SEED';
    # The exact form of the object name list and field list depend on whether
    # we're doing full roles.
    my ($objectNameList, $fieldList);
    if ($roleForm eq 'full') {
        $objectNameList = "Contains HasRole";
        $fieldList = [qw(from-link to-link HasRole(to-link))];
    } else {
        $objectNameList = "Contains";
        $fieldList = [qw(from-link to-link)];
    }
    # Loop through the subsystems. Note we ask Sapling to normalize the
    # subsystem names before we use them.
    for my $subsys (map { $sap->SubsystemID($_) } @$subsystems) {
        # The data for this subsystem will be put in here. It could be
        # a hash or a list.
        my $subsysData;
        # Compute the molecular machine role pattern for this subsystem and
        # genome ID.
        my $pattern = join(":", ERDB::DigestKey($subsys), @patternTail);
        Trace("Machine role pattern is \"$pattern\".") if T(3);
        # Get all the features for the genome in this subsystem. Each fidTuple
        # returned consists of the machine role ID, a feature ID, and a possible
        # role ID. The last piece of the machine role ID is the role abbreviation.
        my @fidTuples = $sap->GetAll($objectNameList, "Contains(from-link) LIKE ?",
                                     [$pattern], $fieldList);
        Trace(scalar(@fidTuples) . " features found.") if T(3);
        # Translate the feature IDs to the desired form. Each feature ID is the
        # second element (index 1) in the tuple.
        for my $fidTuple (@fidTuples) {
            $fidTuple->[1] = $sap->Alias($fidTuple->[1], $source);
        }
        # If we're not sorting by role, then simply store the list.
        if ($roleForm eq 'none') {
            $subsysData = [ map { $_->[1] } @fidTuples ];
            # Convert to a string if we're grouped.
            if ($grouped) {
                $subsysData = join(", ", @$subsysData);
            }

        } else {
            # Here we're going to create a hash of lists, keyed by role. The
            # role ID is determined by the role format.
            $subsysData = {};
            # Loop through the fid tuples.
            for my $fidTuple (@fidTuples) {
                # Get the pieces of data we need from this tuple.
                my ($machineRole, $fid, $role) = @$fidTuple;
                if (! defined $role) {
                    # Here we need to get the role abbreviation from the
                    # machine role ID.
                    (undef, undef, undef, $role) = split /:/, $machineRole;
                }
                # Put this feature in the role hash.
                push @{$subsysData->{$role}}, $fid;
            }
            # Convert the sublists to strings if we're grouped.
            if ($grouped) {
                for my $role (keys %$subsysData) {
                    my $subList = $subsysData->{$role};
                    $subsysData->{$role} = join(", ", @$subList);
                }
            }
        }
        # Store this subsystem's data in the return value.
        $retVal->{$subsys} = $subsysData;
    }
    # Return the result.
    return $retVal;
}

=head3 ids_to_publications

    my $featureHash =       $sapObject->ids_to_publications({
                                -ids => [$id1, $id2, ...],
                                -source => 'UniProt'
    });

Return the PUBMED ID and title of each publication relevant to the specified
feature IDs.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of feature IDs. Normally, these are FIG feature IDs
(e.g. C<fig|100226.1.peg.3361>, C<fig|360108.3.peg.1041>), but other
ID types are permissible if the C<source> parameter is overridden.

=item -source (optional)

Database source of the IDs specified-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
Use C<mixed> to allow mixed ID types (though this may cause problems when the same
ID has different meanings in different databases). Use C<prefixed> to allow IDs with
prefixing indicating the ID type (e.g. C<uni|P00934> for a UniProt ID, C<gi|135813> for
an NCBI identifier, and so forth). The default is C<SEED>.

=item -genome (optional)

ID of a specific genome. If specified, results will only be returned for genes in the
specified genome. The default is to return results for genes in all genomes.

=back

=item RETURN

Returns a reference to a hash mapping feature IDs to lists of 2-tuples. Each
2-tuple consists of a PUBMED ID followed by a publication title.

    $featureHash = { $id1 => [[$pub1a, $title1a], [$pub1b, $title1b], ...],
                     $id2 => [[$pub2a, $title2a], [$pub2b, $title2b], ...],
                     ... };

=back

=cut

sub ids_to_publications {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the attribute database.
    require CustomAttributes;
    my $ca = CustomAttributes->new();
    # Declare the return variable.
    my $retVal = {};
    # Get the list of feature IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Create the feature filter for IDs of this type.
    my ($objects, $filter, @parms) = $sap->ComputeFeatureFilter($args->{-source},
                                                                $args->{-genome});
    # Loop through the identifiers.
    for my $id (@$ids) {
        # The output tuples for this identifier will be put in here.
        my @tuples;
        # Get the proteins for this identifier.
        my @prots = $sap->GetFlat("$objects Produces", $filter, [@parms, $id], 'Produces(to-link)');
        # Get the relevant attributes from the attribute database.
        my @protIDs = map { "Protein:$_" } @prots;
        my @attrs = $ca->GetAttributes(\@protIDs, 'evidence_code', 'dlit%');
        # Convert the attributes into a list of pubmeds for these proteins.
        for my $tuple (@attrs) {
            # Extract the pubmed ID.
            if ($tuple->[2] =~ /dlit\((\d+)/) {
                my $pubmed = $1;
                # Get the publication's hyperlink.
                my ($hyperlink) = $sap->GetEntityValues(Publication => $pubmed, ['citation']);
                # Extract the title.
                my $title = (defined $hyperlink ? $hyperlink->text : "<unknown>");
                # Output the publication ID and title.
                push @tuples, [$pubmed, $title];
            }
        }
        # Store this ID's results.
        $retVal->{$id} = \@tuples;
    }
    # Return the result.
    return $retVal;
}



=head3 ids_to_subsystems

    my $featureHash =       $sapObject->ids_to_subsystems({
                                -ids => [$id1, $id2, ...],
                                -usable => 0,
                                -exclude => ['cluster-based', 'private', ...],
                                -source => 'RefSeq',
                                -subsOnly => 1
                            });

Return the subsystem and role for each feature in the incoming list. A
feature may have multiple roles in a subsystem and may belong to multiple
subsystems, so the role/subsystem information is returned in the form of
a list of ordered pairs for each feature.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of feature IDs. Normally, these are FIG feature IDs
(e.g. C<fig|100226.1.peg.3361>, C<fig|360108.3.peg.1041>), but other
ID types are permissible if the C<source> parameter is overridden.

=item -usable (optional)

If TRUE, then only results from usable subsystems will be included. If FALSE,
then results from all subsystems will be included. The default is TRUE.

=item -exclude (optional)

Reference to a list of special subsystem types that should be excluded from the
result list. The permissible types are C<cluster-based> and C<experimental>.
Normally cluster-based subsystems are included, but experimental subsystems
are only included if the C<-usable> option is turned off.

=item -source (optional)

Database source of the IDs specified-- C<SEED> for FIG IDs, C<GENE> for standard
gene identifiers, or C<LocusTag> for locus tags. In addition, you may specify
C<RefSeq>, C<CMR>, C<NCBI>, C<Trembl>, or C<UniProt> for IDs from those databases.
Use C<mixed> to allow mixed ID types (though this may cause problems when the same
ID has different meanings in different databases). Use C<prefixed> to allow IDs with
prefixing indicating the ID type (e.g. C<uni|P00934> for a UniProt ID, C<gi|135813> for
an NCBI identifier, and so forth). The default is C<SEED>.

=item -genome (optional)

ID of a specific genome. If specified, results will only be returned for genes in the
specified genome. The default is to return results for genes in all genomes.

=item -subsOnly (optional)

If TRUE, instead of a list of (role, subsystem) 2-tuples, each feature ID will be
mapped to a simple list of subsystem names. The default is FALSE.

=back

=item RETURN

Returns a reference to a hash mapping feature IDs to lists of 2-tuples. Each
2-tuple consists of a role name followed by a subsystem name. If a feature is
not in a subsystem, it will not be present in the return hash.

=over 8

=item Normal Output

    $featureHash = { $id1 => [[$role1a, $sub1a], [$role1b, $sub1b], ...],
                     $id2 => [[$role2a, $sub2a], [$role2b, $sub2b], ...],
                     ... };

=item Output if -subsOnly = 1

    $featureHash = { $id1 => [$sub1a, $sub1b, ...],
                     $id2 => [$sub2a, $sub2b, ...],
                     ... };

=back

=back

=cut

sub ids_to_subsystems {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Get the desired ID type.
    my $source = $args->{-source};
    # Build the filter string, object name list, and parameter value list prefix.
    my ($objects, $filter, @parms) = $sap->ComputeFeatureFilter($source,
                                                                $args->{-genome});
    ServerThing::AddSubsystemFilter(\$filter, $args);
    # Loop through the features.
    for my $id (@$ids) {
        Trace("Looking for subsystems of feature $id.") if T(3);
        # Get this feature's subsystem information.
        my @rows = $sap->GetAll("$objects IsContainedIn MachineRole IsRoleFor Implements Variant IsDescribedBy Subsystem AND MachineRole HasRole",
                                 $filter, [@parms, $id], [qw(HasRole(to-link)
                                                             IsDescribedBy(to-link))]);
        # Only proceed if results were found.
        if (@rows) {
            # Determine whether we want the roles and the names or just the names. Whatever
            # we do need will be put in here.
            my $subData;
            if ($args->{-subsOnly}) {
                # Here we want only the subsystem names. We use a hash to remove
                # duplicates.
                my %subNames = map { $_->[1] => 1 } @rows;
                $subData = [ sort keys %subNames ];
            } else {
                # Here we want everything.
                $subData = \@rows;
            }
            # Store it in the result hash.
            $retVal->{$id} = $subData;
        }
    }
    # Return the result.
    return $retVal;
}


=head3 is_in_subsystem

    my $featureHash =       $sapObject->is_in_subsystem({
                                -ids => [$fid1, $fid2, ...],
                                -usable => 0,
                                -exclude => [$type1, $type2, ...]
                            });

Return the subsystem and role for each specified feature.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the FIG feature IDs for the features of interest.

=item -usable (optional)

If TRUE, then only results from usable subsystems will be included. If FALSE,
then results from all subsystems will be included. The default is TRUE.

=item -exclude (optional)

Reference to a list of special subsystem types that should be excluded from the
result list. The permissible types are C<cluster-based> and C<experimental>.
Normally cluster-based subsystems are included, but experimental subsystems
are only included if the C<-usable> option is turned off.

=back

For backward compatibility, the parameter may also be a reference to a list
of FIG feature IDs.

=item RETURN

Returns a reference to a hash that maps each incoming feature ID to a list of
2-tuples, each 2-tuple consisting of (0) the ID of a subsystem containing the
feature and (1) the feature's role in that subsystem. If an incoming feature is
not in any subsystem, its ID will be mapped to an empty list.

    $featureHash = { $fid1 => [[$sub1a, $role1a], [$sub1b, $role1b], ...],
                     $fid2 => [[$sub2a, $role2a], [$sub2b, $role2b[, ...],
                     ... };

In backward-compatible mode, returns a reference to a list of 3-tuples, each
3-tuple consisting of (0) a subsystem ID, (1) a role ID, and (2) the ID of a
feature from the input list.

=back

=cut

sub is_in_subsystem {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # This will be set to TRUE if we are in backward-compatible mode.
    my $backwardMode = 0;
    # Convert a list to a hash.
    if (ref $args ne 'HASH') {
        $args = { -ids => $args };
        $backwardMode = 1;
    }
    # Create the filter clause. It contains at least a feature filter.
    my $filter = 'Feature(id) = ?';
    # There may also be subsystem filtering.
    ServerThing::AddSubsystemFilter(\$filter, $args);
    # Declare the return variable.
    my $retVal = {};
    # Get the fig IDs from the parameters.
    my $ids = ServerThing::GetIdList(-ids => $args);
    foreach my $fid (@$ids) {
        my @resultRows = $sapling->GetAll("Feature IsContainedIn MachineRole HasRole Role AND " .
                                          "MachineRole IsRoleFor MolecularMachine Implements Variant IsDescribedBy Subsystem",
                                          $filter, [$fid], [qw(Subsystem(id) Role(id))]);
        $retVal->{$fid} = \@resultRows;
    }
    # If we're in backward-compatible mode, convert the return value to a list.
    if ($backwardMode) {
        my @list;
        for my $fid (@$ids) {
            push @list, map { [@$_, $fid] } @{$retVal->{$fid}};
        }
        $retVal = \@list;
    }
    # Return the result.
    return $retVal;
}

=head3 is_in_subsystem_with

    my $featureHash =       $sapObject->is_in_subsystem_with({
                                -ids => [$fid1, $fid2, ...],
                                -usable => 0,
                                -exclude => [$type1, $type2, ...]
                            });

For each incoming feature, returns a list of the features in the same genome that
are part of the same subsystem. For each other feature returned, its role,
functional assignment, subsystem variant, and subsystem ID will be returned as
well.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the FIG feature IDs for the features of interest.

=item -usable (optional)

If TRUE, then only results from usable subsystems will be included. If FALSE,
then results from all subsystems will be included. The default is TRUE.

=item -exclude (optional)

Reference to a list of special subsystem types that should be excluded from the
result list. The permissible types are C<cluster-based> and C<experimental>.
Normally cluster-based subsystems are included, but experimental subsystems
are only included if the C<-usable> option is turned off.

=back

For backward compatibility, the parameter may also be a reference to a list
of FIG feature IDs.

=item RETURN

Returns a reference to a hash that maps each incoming feature ID to a list
of 5-tuples relating to features in the same subsystem. Each 5-tuple contains
(0) a subsystem ID, (1) a variant ID, (2) the related feature ID, (3) the
related feature's functional assignment, and (4) the related feature's role
in the subsystem.

    $featureHash = { $fid1 => [[$sub1a, $variant1a, $fid1a, $function1a, $role1a],
                               [$sub1b, $variant1b, $fid1b, $function1b, $role1b], ...],
                     $fid2 => [[$sub2a, $variant2a, $fid2a, $function2a, $role2a],
                               [$sub2b, $variant2b, $fid2b, $function2b, $role2b], ...],
                    ... };

In backward-compatibility mode, returns a reference to a list of lists. Each
sub-list contains 6-tuples relating to a single incoming feature ID. Each
6-tuple consists of a subsystem ID, a variant ID, the incoming feature ID, the
other feature ID, the other feature's functional assignment, and the other
feature's role in the subsystem.

=back

=cut

sub is_in_subsystem_with {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # Declare the return variable.
    my $retVal;
    # This will be set to TRUE if we are in backward-compatible mode.
    my $backwardMode = 0;
    # Convert a list to a hash.
    if (ref $args ne 'HASH') {
        $args = { -ids => $args };
        $backwardMode = 1;
    }
    # Create the filter clause. It contains at least a feature filter.
    my $filter = 'Feature(id) = ?';
    # Add any required subsystem filtering.
    ServerThing::AddSubsystemFilter(\$filter, $args);
    # Get the fig IDs from the parameters.
    my $ids = ServerThing::GetIdList(-ids => $args);
    foreach my $fid (@$ids) {
        my @resultRows = $sapling->GetAll("Feature IsContainedIn MachineRole IsRoleFor MolecularMachine Implements Variant IsDescribedBy Subsystem AND MolecularMachine IsMachineOf MachineRole2 Contains Feature2 AND MachineRole2 HasRole Role",
                                          $filter, [$fid],
                                          [qw(Subsystem(id) Variant(code)
                                              Feature2(id) Feature2(function)
                                              Role(id))]);
        $retVal->{$fid} = \@resultRows;
    }
    # If this is backward-compatability mode, convert the result to a list.
    if ($backwardMode) {
        my @outList;
        for my $fid (@$ids) {
            my $fidList = $retVal->{$fid};
            if (! defined $fidList) {
                push @outList, [];
            } else {
                # Because the incoming feature ID is no longer available as the
                # hash key, we need to put it back into the output tuples. It goes
                # in the third position (offset 2).
                for my $fidTuple (@$fidList) {
                    splice @$fidTuple, 2, 0, $fid;
                }
                push @outList, $fidList;
            }
        }
        $retVal = \@outList;
    }
    # Return the result.
    return $retVal;
}

=head3 pegs_implementing_roles

    my $roleHash =          $sapObject->pegs_implementing_roles({
                                -subsystem => $subsysID,
                                -roles => [$role1, $role2, ...]
                            });

Given a subsystem and a list of roles, return a list of the subsystem's
features for each role.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 4

=item -subsystem

ID of a subsystem.

=item -roles

Reference to a list of roles.

=back

For backward compatibility, the parameter can also be a reference to a 2-tuple
consisting of (0) a subsystem ID and (1) a reference to a list of roles.

=item RETURN

Returns a hash that maps each role ID to a list of the IDs for the features that
perform the role in that subsystem.

    $roleHash = { $role1 => [$fid1a, $fid1b, ...],
                  $role2 => [$fid2a, $fid2b, ...],
                  ... };

In backward-compatibility mode, returns a list of 2-tuples. Each tuple consists
of a role and a reference to a list of the features in that role.

=back

=cut

sub pegs_implementing_roles {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # Get the sapling subsystem object.
    require SaplingSubsys;
    # Declare the return variable.
    my $retVal = {};
    # Check for backward-compatibility mode.
    my $backwardMode = 0;
    if (ref $args ne 'HASH') {
        $args = { -subsystem => $args->[0], -roles => $args->[1] };
        $backwardMode = 1;
    }
    # Get the subsystem ID.
    my $subsystem = $args->{-subsystem};
    # Get the list of roles.
    my $roles = ServerThing::GetIdList(-roles => $args);
    # If there is no subsystem ID, it's an error.
    if (! defined $subsystem) {
        Confess("Subsystem ID not specified.");
    } else {
        # Normalize the subsystem ID.
        my $subsystemID = $sapling->SubsystemID($subsystem);
        # Get a sapling subsystem object.
        my $ss = SaplingSubsys->new($subsystemID, $sapling);
        # Only proceed if we found one.
        if (defined $ss) {
            # Loop through the roles, acquiring features.
            foreach my $role (@$roles) {
                $retVal->{$role} = [$ss->pegs_for_role($role)];
            }
        }
    }
    # In backward-compatible mode, we must convert the output to a list.
    if ($backwardMode) {
        my @outList = map { [$_, $retVal->{$_} ] } @$roles;
        $retVal = \@outList;
    }
    # Return the result.
    return $retVal;
}

=head3 pegs_in_subsystems

    my $subsysHash =        $sapObject->pegs_in_subsystems({
                                -genomes => [$genome1, $genome2, ...],
                                -subsystems => [$sub1, $sub2, ...]
                            });

This method takes a list of genomes and a list of subsystems and returns
a list of the roles represented in each genome/subsystem pair.

=over 4

=item parameter

Reference to a hash of parameter values with the following possible keys.

=over 8

=item -genomes

Reference to a list of genome IDs.

=item -subsystems

Reference to a list of subsystem IDs.

=back

For backward compatibility, the parameter may also be a reference to a 2-tuple,
the first element of which is a list of genome IDs and the second of which is a
list of subsystem IDs.

=item RETURN

Returns a reference to a hash of hashes. The main hash is keyed by subsystem ID.
Each subsystem's hash is keyed by role ID and maps the role to a list of
the feature IDs for that role in the subsystem that belong to the specified
genomes.

    $subsysHash = { $sub1 => { $role1A => [$fid1Ax, $fid1Ay, ...],
                               $role1B => [$fid1Bx, $fid1By, ...],
                               ... },
                    $sub2 => { $role2A => [$fid2Ax, $fid2Ay, ...],
                               $role2B => [$fid2Bx, $fid2By, ...],
                               ... },
                    ... };

In backward-compatibility mode, returns a list of 2-tuples. Each tuple consists
of a subsystem ID and a second 2-tuple that contains a role ID and a reference
to a list of the feature IDs for that role that belong to the specified genomes.

=back

=cut

sub pegs_in_subsystems {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # Get access to the sapling subsystem object.
    require SaplingSubsys;
    # Declare the return variable.
    my $retVal = {};
    # Check for backward-compatibility mode.
    my $backwardMode = 0;
    if (ref $args ne 'HASH') {
        $args = { -genomes => $args->[0], -subsystems => $args->[1] };
        $backwardMode = 1;
    }
    # Get the list of genome IDs.
    my $genomes = ServerThing::GetIdList(-genomes => $args);
    # Get the list of subsystem IDs.
    my $subs = ServerThing::GetIdList(-subsystems => $args);
    # Loop through the subsystems.
    for my $sub (@{$subs}) {
        # Normalize the subsystem ID.
        my $subID = $sapling->SubsystemID($sub);
        # Get the subsystem spreadsheet in memory.
        my $ss = SaplingSubsys->new($subID, $sapling);
        # Only proceed if we found it.
        if (defined $ss) {
            # We'll build the subsystem's hash in here.
            my $subHash = {};
            # Loop through the genomes, assigning features to the roles.
            foreach my $g (@{$genomes}) {
                # Get role/featureList pairs for this genome.
                my @roleTuples = $ss->get_roles_for_genome($g, 1);
                # Loop through the pairs.
                foreach my $roleTuple (@roleTuples) {
                    # Extract the role ID and the feature list.
                    my ($role, $features) = @$roleTuple;
                    # Attach the features to the role.
                    push @{$subHash->{$role}}, @$features;
                }
            }
            # Attach this hash to this subsystem.
            $retVal->{$sub} = $subHash;
        }
    }
    # In backward-compatible mode, we have to convert the hashes to lists.
    if ($backwardMode) {
        # We'll build the output list in here.
        my @outList;
        # Loop through the subsystems in input order.
        for my $ss (@$subs) {
            my $subHash = $retVal->{$ss};
            if (defined $subHash) {
                # Now we convert the role -> feature map to a list of
                # [sub, [role, feature]] nested pairs.
                for my $role (keys %$subHash) {
                    push @outList, [$ss, [$role, $subHash->{$role}]];
                }
            }
        }
        # Store the output list as the result.
        $retVal = \@outList;
    }
    # Return the result.
    return $retVal;
}

# Synonym for "pegs_in_subsystems" provided for backward compatibility.
sub pegs_in_subsystem {
    return pegs_in_subsystems(@_);
}

=head3 pegs_in_variants

    my $subsysHash =        $sapObject->pegs_in_variants({
                                -genomes => [$genomeA, $genomeB, ...],
                                -subsystems => [$sub1, $sub2, ...]
                            });

This method takes a list of genomes and a list of subsystems and returns
a list of the pegs represented in each genome/subsystem pair.

The main difference between this method and L</pegs_in_subsystems> is in
the organization of the output, which is more like a subsystem spreadsheet.

=over 4

=item parameter

Reference to a hash of parameter values with the following possible keys.

=over 8

=item -genomes (optional)

Reference to a list of genome IDs. If the list is omitted, all genomes will be
included in the output (which will be rather large in most cases).

=item -subsystems

Reference to a list of subsystem IDs.

=back

=item RETURN

Returns a reference to a hash mapping subsystem IDs to sub-hashes. Each sub-hash
is keyed by genome ID and maps the genome ID to a list containing the variant code
and one or more n-tuples, each n-tuple containing a role ID followed by a list of
the genes in the genome having that role in the subsystem.

    $subsysHash = { $sub1 => { $genomeA => [$vc1A,
                                            [$role1Ax, $fid1Ax1, $fid1Ax2, ...],
                                            [$role1Ay, $fid1Ay1, $fid1Ay2, ...],
                                            ...],
                               $genomeB => [$vc1B,
                                            [$role1Bx, $fid1Bx1, $fid1Bx2, ...],
                                            [$role1By, $fid1By1, $fid1By2, ...],
                                            ...],
                               ... },
                    $sub2 => { $genomeA => [$vc2A,
                                            [$role2Ax, $fid2Ax1, $fid2Ax2, ...],
                                            [$role2Ay, $fid2Ay1, $fid2Ay2, ...],
                                            ...],
                               $genomeB => [$vc2B,
                                            [$role2Bx, $fid2Bx1, $fid2Bx2, ...],
                                            [$role2By, $fid2By1, $fid2By2, ...],
                                            ...],
                               ... },
                    ... };

Note that in some cases the genome ID will include a region string. This happens
when the subsystem has multiple occurrences in the genome.

=back

=cut

sub pegs_in_variants {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the genome IDs.
    my $genomes = ServerThing::GetIdList(-genomes => $args, 1);
    # Form them into a hash.
    my %genomeH = map { $_ => 1 } @$genomes;
    # If no genomes were specified, we allow all of them.
    my $allGenomes = (! @$genomes);
    # Get the subsystem IDs.
    my $subs = ServerThing::GetIdList(-subsystems => $args);
    # Loop through the subsystems.
    for my $sub (@$subs) {
        # This will be the hash for this subsystem.
        my %subHash;
        # Get this subsystem's molecular machines.
        my @machines = $sap->GetAll("Describes Variant IsImplementedBy MolecularMachine IsUsedBy",
                                    'Describes(from-link) = ?', [$sub], [qw(Variant(code)
                                    MolecularMachine(id) MolecularMachine(region)
                                    IsUsedBy(to-link))]);
        # Loop through the machines, pausing on the ones related to our genomes.
        for my $machine (@machines) {
            # Get this machine's data.
            my ($variant, $machineID, $region, $genome) = @$machine;
            # Only proceed if it's a genome of interest.
            if ($allGenomes || $genomeH{$genome}) {
                # Get the features and roles for this machine.
                my @ssData = $sap->GetAll("IsMachineOf MachineRole HasRole AND MachineRole Contains",
                                          'IsMachineOf(from-link) = ?', [$machineID],
                                          [qw(HasRole(to-link) Contains(to-link))]);
                # Only proceed if we found some.
                if (@ssData) {
                    # Create the full genome ID using the region.
                    my $genomeID = $genome;
                    if ($region) {
                        $genomeID .= ":$region";
                    }
                    # The cells will be accumulated in here.
                    my @row = $variant;
                    # Get the first role.
                    my $role = $ssData[0][0];
                    # This will accumulate the current cell.
                    my @cell;
                    for my $ssTuple (@ssData) {
                        my ($newRole, $fid) = @$ssTuple;
                        # If this is a new role, output the old cell.
                        if ($newRole ne $role) {
                            push @row, [$role, @cell];
                            $role = $newRole;
                            @cell = $fid;
                        } else {
                            # It's the old role, so add the gene to the cell.
                            push @cell, $fid;
                        }
                    }
                    # Store this genome's data in the subsystem's hash.
                    $subHash{$genomeID} = [@row, [$role, @cell]];
                }
            }
        }
        # Store this subsystem's data in the return hash.
        $retVal->{$sub} = \%subHash;
    }
    # Return the result.
    return $retVal;
}

=head3 roles_exist_in_subsystem

    my $rolesHash =         $sapObject->roles_exist_in_subsystem({
                                -subsystem => $sub1,
                                -roles => [$role1, $role2, ...]
                            });

Indicate which roles in a given list belong to a specified subsystem.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -subsystem

The name of the subsystem of interest.

=item -roles

A reference to a list of role IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming role ID to C<1> if it
exists in the specified subsystem and C<0> otherwise.

    $roleHash = { $role1 => $flag1, $role2 => $flag2, ... };

=back

=cut

sub roles_exist_in_subsystem {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Declare the return variable.
    my $retVal = {};
    # Get the list of roles.
    my $roles = ServerThing::GetIdList(-roles => $args);
    # Get the subsystem name.
    my $subsystem = $args->{-subsystem};
    Confess("No -subsystem specified.") if ! defined $subsystem;
    # Normalize the subsystem name.
    $subsystem = $sap->SubsystemID($subsystem);
    # Create a hash of all the roles in the subsystem.
    my %roleHash = map { $_ => 1 } $sap->GetFlat("Includes",
                                                 "Includes(from-link) = ?",
                                                 [$subsystem],
                                                 'Includes(to-link)');
    # Loop through the incoming roles, comparing them against the
    # hash.
    for my $role (@$roles) {
        $retVal->{$role} = ($roleHash{$role} ? 1 : 0);
    }
    # Return the result.
    return $retVal;
}

=head3 roles_to_subsystems

    my $roleHash =              $sapObject->({
                                    -roles => [$role1, $role2, ...],
                                    -usable => 0
                                });

Return the subsystems containing each specified role.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -roles

Reference to a list of role names.

=item -usable (optional)

If TRUE, only usable subsystems will be returned. If FALSE, all subsystems
will be returned. The defult is TRUE.

=back

=item RETURN

Returns a reference to a hash mapping each incoming role to a list of
the names of subsystems containing that role.

    $roleHash = { $role1 => [$sub1a, $sub1b, ...],
                  $role2 => [$sub2a, $sub2b, ...],
                  ... };

=back

=cut

sub roles_to_subsystems {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the Sapling database.
    my $sap = $self->{db};
    # Declare the return hash.
    my $retVal = {};
    # Get the incoming roles.
    my $roles = ServerThing::GetIdList(-roles => $args);
    # Get the usability flag.
    my $usable = $args->{-usable} || 0;
    # Create the filter clause.
    my $filter = "IsIncludedIn(from-link) = ?";
    if ($usable) {
        $filter .= " AND Subsystem(usable) = 1";
    }
    # Loop through the roles.
    for my $role (@$roles) {
        # Get the subsystems for this role.
        my (@subs) = $sap->GetFlat("IsIncludedIn Subsystem", $filter,
            [$role], 'Subsystem(id)');
        # Store them in the return hash.
        $retVal->{$role} = \@subs;
    }
    # Return the result hash.
    return $retVal;
}

=head3 rows_of_subsystem

    my $subHash =               $sapObject->({
                                    -subs => [$sub1, $sub2, ...],
                                    -genomes => [$genomeA, $genomeB, ...],
                                    ...
                                });

Return the subsystem row for each subsystem/genome pair. A row in this
case consists of a reference to a hash mapping role names to a list of
the FIG feature IDs for the features in the genome performing that
role.

In the Sapling database, a subsystem row is represented by the
B<MolecularMachine> entity. The strategy of this method is therefore
to find the molecular machine for each subsystem/genome pair, and
then use its ID to get the roles and features.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -subs

Reference to a list of subsystem IDs.

=item -genomes

Reference to a list of genome IDs.

=back

=item RETURN

Returns a reference to a hash mapping each incoming subsystem ID to
a sub-hash keyed by genome ID. In the sub-hash, each genome ID will
map to a sub-sub-hash that maps role names to lists of feature IDs.

    $subHash = { $sub1 => { $genomeA => { $role1Aa => [$fid1Aax, $fid1Aay, ... ],
                                          $role1Ab => [$fid1Abx, $fid1Aby, ... ],
                                          ... },
                            $genomeB => { $role1Ba => [$fid1Bax, $fid1Bay, ... ],
                                          $role1Bb => [$fid1Bbx, $fid1Bby, ... ],
                                          ... },
                            ... },
                 $sub2 => { $genomeA => { $role2Aa => [$fid2Aax, $fid2Aay, ... ],
                                          $role2Ab => [$fid2Abx, $fid2Aby, ... ],
                                          ... },
                            $genomeB => { $role2Ba => [$fid2Bax, $fid2Bay, ... ],
                                          $role2Bb => [$fid2Bbx, $fid2Bby, ... ],
                                          ... },
                            ... },
                 ... };

=back

=cut

sub rows_of_subsystems {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of subsystems.
    my $subs = ServerThing::GetIdList(-subs => $args);
    # Get the list of genomes.
    my $genomes = ServerThing::GetIdList(-genomes => $args);
    # Declare the return hash.
    my $retVal = {};
    # Loop through the subsystems.
    for my $subID (@$subs) {
        # Normalize the subsystem ID.
        my $sub = $sap->SubsystemID($subID);
        # This hash will map genomes to rows.
        my %rows;
        # Loop through the genomes.
        for my $genome (@$genomes) {
            # This hash will map roles to feature lists.
            my %cells;
            # Get the molecular machine ID for this subsystem and genome.
            my ($rowID) = $sap->GetFlat("Describes IsImplementedBy MolecularMachine IsUsedBy",
                            'Describes(from-link) = ? AND IsUsedBy(to-link) = ?',
                            [$sub,$genome], 'MolecularMachine(id)');
            # If we found it (there can be at most one), get the role/feature
            # pairs.
            if ($rowID) {
                my @pairs = $sap->GetAll("IsMachineOf MachineRole HasRole AND MachineRole Contains",
                                'IsMachineOf(from-link) = ?', [$rowID],
                                [qw(HasRole(to-link) Contains(to-link))]);
                # Form them into a hash of roles to feature lists.
                for my $pair (@pairs) {
                    push @{$cells{$pair->[0]}}, $pair->[1];
                }
            }
            # Store this row in the hash of rows for this subsystem.
            $rows{$genome} = \%cells;
        }
        # Store the rows for this subsystem in the return hash.
        $retVal->{$subID} = \%rows;
    }
    # Return the result hash.
    return $retVal;
}


=head3 subsystem_data

    my $subsysHash =        $sapObject->subsystem_data({
                                -ids => [$sub1, $sub2, ...],
                                -field => 'version'
                            });

For each incoming subsystem ID, return the specified data field. This
method can be used to find the curator, description, or version of the
specified subsystems.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of subsystem IDs.

=item -field (optional)

Name of the desired data field-- C<curator> to retrieve the name of each
subsystem's curator, C<version> to get the subsystem's version number,
or C<description> to get the subsystem's description, or C<notes> to
get the subsystem's notes. The default is C<description>.

=back

=item RETURN

Returns a hash mapping each incoming subsystem ID to the associated data
value.

    $subsysHash = { $sub1 => $value1, $sub2 => $value2, ... };

=back

=cut

use constant VALID_SUB_DATA_FIELDS => { curator     => 1,
                                        description => 1,
                                        version     => 1,
                                        notes       => 1
                                      };

sub subsystem_data {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the name of the desired field and insure it's valid.
    my $field = $args->{-field} || 'description';
    if (! VALID_SUB_DATA_FIELDS->{$field}) {
        Confess("Invalid subsystem field \"$field\" specified.");
    } else {
        # Get the list of subsystem IDs.
        my $ids = ServerThing::GetIdList(-ids => $args);
        # Loop through the subsystems, retrieving the desired data field.
        for my $id (@$ids) {
            # Normalize the ID.
            my $realID = $sap->SubsystemID($id);
            # Get the desired value.
            my ($value) = $sap->GetEntityValues(Subsystem => $realID, [$field]);
            # If we found something, put it in the return hash.
            if (defined $value) {
                $retVal->{$id} = $value;
            }
        }
    }
    # Return the result.
    return $retVal;
}

=head3 subsystem_genomes

    my $subHash =           $sapObject->subsystem_genomes({
                                -ids => [$sub1, $sub2, ...],
                                -all => 1
                            });

For each subsystem, return the genomes that participate in it and their
associated variant codes.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the names of the subsystems whose genome information is
desired.

=item -all (optional)

If TRUE, then all genomes associated with the subsystem will be listed. The
default is FALSE, meaning that only genomes that completely implement the
subsystem will be listed.

=back

=item RETURN

Returns a reference to a hash that maps each subsystem ID to a sub-hash. Each
sub-hash in turn maps the ID of each subsystem that participates in the subsystem
to its variant code.

    $subHash = { $sub1 => { $genome1a => $code1a, $genome1b => $code1b, ...},
                 $sub2 => { $genome2a => $code2a, $genome2b => $code2b, ...},
                 ... };

=back

=cut

sub subsystem_genomes {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of subsystem IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Compute the filter clause for the variants.
    my $filter = 'Describes(from-link) = ?';
    if (! $args->{-all}) {
        $filter .= "AND Variant(type) = 'normal'";
    }
    # Declare the return variable.
    my $retVal = {};
    # Loop through the subsystem IDs.
    for my $id (@$ids) {
        # Get the genome and variant information.
        my %genomes = map { $_->[0] => $_->[1] }
                        $sap->GetAll("Describes Variant IsImplementedBy IsUsedBy",
                        $filter, [$id], [qw(IsUsedBy(to-link) Variant(code))]);
        # Store it in the result hash.
        $retVal->{$id} = \%genomes;
    }
    # Return the result.
    return $retVal;
}


=head3 subsystem_names

    my $nameList =          $sapObject->subsystem_names({
                                -usable => 0,
                                -exclude => ['cluster-based', ...]
                            });

Return a list of all subsystems in the database.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -usable (optional)

If TRUE, then only results from usable subsystems will be included. If FALSE,
then results from all subsystems will be included. The default is TRUE.

=item -exclude (optional)

Reference to a list of special subsystem types that should be excluded from the
result list. The permissible types are C<cluster-based> and C<experimental>.
Normally cluster-based subsystems are included, but experimental subsystems
are only included if the C<-usable> option is turned off.

=back

=item RETURN

Returns a reference to a list of subsystem names.

    $nameList = [$sub1, $sub2, ...];

=back

=cut

sub subsystem_names {
    # Get the parameters.
    my ($self, $args) = @_;
    # Form the filter string.
    my $filter = "";
    ServerThing::AddSubsystemFilter(\$filter, $args);
    # Ask for the subsystem IDs.
    my $retVal = [ $self->{db}->GetFlat("Subsystem", $filter, [], 'id') ];
    # Return the result.
    return $retVal;
}

=head3 subsystem_roles

    my $subHash =           $sapObject->subsystem_roles({
                                -ids => [$sub1, $sub2, ...],
                                -aux => 1
    });

Return the list of roles for each subsystem, in order.

=over 4

=item parameter

Reference to a hash of parameters with the following possible keys.

=over 8

=item -ids

Reference to a list of subsystem IDs.

=item -aux (optional)

If TRUE, auxiliary roles will be included. The default is FALSE,
which excludes auxiliary roles.

=item -abbr (optional)

If TRUE, then the role abbreviations will be included in the results. In
this case, each subsystem name will be mapped to a list of 2-tuples, with
each 2-tuple consisting of (0) the role name and (1) the role abbreviation.
The default is FALSE (normal output).

=back

=item RETURN

Return a hash mapping each subsystem ID to a list of roles (normal) or
a list of role/abbreviation pairs (extended output).

=over 8

=item Output if -abbr is FALSE

    $subHash = { $sub1 => [$role1a, $role1b, ...],
                 $sub2 => [$role2a, $role2b, ...],
                 ... };

=item Output if -abbr is TRUE

    $subHash = { $sub1 => [[$role1a, $abbr1a],
                           [$role1b, $abbr1b], ...],
                 $sub2 => [[$role2a, $abbr2a],
                           [$role2b, $abbr2b], ...],
                 ... };

=back

=back

=cut

sub subsystem_roles {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sap = $self->{db};
    # Get the list of subsystem IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Compute the filter.
    my $filter = 'Subsystem(id) = ?';
    if (! $args->{-aux}) {
        $filter .= ' AND Includes(auxiliary) = 0';
    }
    $filter .= ' ORDER BY Includes(sequence)';
    # Declare the return variable.
    my $retVal = {};
    # Loop through the IDs.
    for my $id (@$ids) {
        # Normalize the subsystem ID.
        my $realID = $sap->SubsystemID($id);
        # We'll put the roles for this subsystem in here.
        my @results;
        # Process according to the type of output desired.
        if ($args->{-abbr}) {
            @results = $sap->GetAll("Subsystem Includes Role", $filter,
                                     [$realID],
                                     'Role(id) Includes(abbreviation)');
        } else {
            @results = $sap->GetFlat("Subsystem Includes Role", $filter,
                                     [$realID], 'Role(id)');
        }
        # Store them in the results.
        $retVal->{$id} = \@results;
    }
    # Return the result hash.
    return $retVal;
}

=head3 subsystem_spreadsheet

    my $subsysHash =        $sapObject->subsystem_spreadsheet({
                                -ids => [$sub1, $sub2, ...]
                            });

This method takes a list of subsystem IDs, and for each one returns a
list of the features in the subsystem. For each feature, it will include
the feature's functional assignment, the subsystem name and variant
(spreadsheet row), and its role (spreadsheet column).

=over 4

=item parameter

Reference to a hash of parameters with the following possible keys.

=over 8

=item -ids

Reference to a list of subsystem IDs.

=back

For backward compatibility, this method can also accept a reference to a list of
subsystem IDs.

=item RETURN

Returns a hash mapping each incoming subsystem ID to a list of 4-tuples. Each
tuple contains (0) a variant ID, (1) a feature ID, (2) the feature's functional
assignment, and (3) the feature's role in the subsystem.

    $subsysHash = { $sub1 => [[$variant1a, $fid1a, $function1a, $role1a],
                              [$variant1b, $fid1b, $function1b, $role1b], ...],
                    $sub2 => [[$variant2a, $fid2a, $function2a, $role2a],
                              [$variant2b, $fid2b, $function2b, $role2b], ...],
                    ... };

In backward-compatability mode, returns a list of 5-tuples. Each tuple contains
(0) a subsystem ID, (1) a variant ID, (2) a feature ID, (3) the feature's
functional assignment, and (4) the feature's role in the subsystem.

=back

=cut

sub subsystem_spreadsheet {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # Declare the return variable.
    my $retVal;
    # Check for the backward-compatible mode.
    my $backwardMode = 0;
    if (ref $args ne 'HASH') {
        $args = { -ids => $args };
        $backwardMode = 1;
    }
    # Get the list of subsystem IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Loop through the subsystem IDs.
    foreach my $subsysName (@$ids) {
        # Normalize the subsystem ID.
        my $subsysID = $sapling->SubsystemID($subsysName);
        # Get the subsystem's spreadsheet data.
        my @resultRows = $sapling->GetAll("Subsystem Describes Variant IsImplementedBy MolecularMachine IsMachineOf MachineRole Contains Feature AND MachineRole HasRole Role Includes Subsystem",
                                          'Subsystem(id) = ? ORDER BY Variant(id), Includes(sequence)',
                                          [$subsysID], [qw(Variant(id)
                                                           Feature(id)
                                                           Feature(function)
                                                           Role(id))]);
        $retVal->{$subsysName} = \@resultRows;
    }
    # In backward-compatible mode, convert the hash to a list.
    if ($backwardMode) {
        # We'll build the list in here.
        my @listForm;
        for my $subsysName (@$ids) {
            # Get this subsystem's spreadsheet and paste in the subsystem ID.
            my $spreadsheet = $retVal->{$subsysName};
            for my $row (@$spreadsheet) {
                unshift @$row, $subsysName;
            }
            # Put it into the output.
            push @listForm, @$spreadsheet;
        }
        # Return the list.
        $retVal = \@listForm;
    }
    # Return the result.
    return $retVal;
}


=head3 subsystem_type

    my $subsysHash =        $sapObject->subsystem_type({
                                -ids => [$sub1, $sub2, ...],
                                -type => 'cluster-based'
                            });

For each incoming subsystem, return TRUE if it has the specified
characteristic, else FALSE.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of subsystem names.

=item -type

Name of the subsystem characteristic of interest. The default is C<usable>. The
possible characteristics are

=over 8

=item cluster-based

A I<cluster-based> subsystem is one in which there is functional-coupling
evidence that genes belong together, but we do not yet know what they do.

=item experimental

An I<experimental> subsystem is designed for investigation and is not yet ready to
be used in comparative analysis and annotation.

=item private

A I<private> subsystem has valid data, but is not considered ready for general
distribution.

=item usable

An unusable subsystem is one that is experimental or is of such low quality that
it can negatively affect analysis. A I<usable> subsystem is one that is not
unusable.

=back

=back

=item RETURN

Returns a hash mapping the incoming subsystem names to TRUE/FALSE flags indicating
the value of the specified characteristic.

    $subsysHash = { $sub1 => $flag1, $sub2 => $flag2, ... };

=back

=cut

use constant VALID_SUBSYSTEM_TYPES => { 'cluster-based' => 1,
                                        experimental    => 1,
                                        private         => 1,
                                        usable          => 1
                                      };

sub subsystem_type {
    # Get the parameters.
    my ($self, $args) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get the Sapling database.
    my $sap = $self->{db};
    # Get the subsystem IDs.
    my $ids = ServerThing::GetIdList(-ids => $args);
    # Get the characteristic of interest. Note the default is "usable".
    my $type = $args->{-type} || 'usable';
    # Insure it's valid.
    if (! VALID_SUBSYSTEM_TYPES->{$type}) {
        Confess("Invalid subsystem type \"$type\" specified.");
    } else {
        # Loop through the subsystem IDs.
        for my $id (@$ids) {
            # Normalize the ID.
            my $realID = $sap->SubsystemID($id);
            # Get the desired flag value.
            my ($flag) = $sap->GetEntityValues(Subsystem => $realID, [$type]);
            # Store it in the return hash.
            $retVal->{$id} = $flag;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 subsystems_for_role

    my $roleHash =          $sapObject->subsystems_for_role({
                                -ids => [$role1, $role2, ...],
                                -usable => 1,
                                -exclude => ['cluster-based', ...]
                            });

For each role, return a list of the subsystems containing that role. The
results can be filtered to include unusable subsystems or exclude subsystems
of certain exotic types.

=over 4

=item parameter

The parameter should be a reference to a hash with the following keys.

=over 8

=item -ids

Reference to a list of the IDs of the roles of interest.

=item -aux (optional)

If TRUE, then subsystems in which the role is auxiliary will be included.
The default is not to include such subsystems.

=item -usable (optional)

If TRUE, then only results from usable subsystems will be included. If FALSE,
then results from all subsystems will be included. The default is TRUE.

=item -exclude (optional)

Reference to a list of special subsystem types that should be excluded from the
result list. The permissible types are C<cluster-based> and C<experimental>.
Normally cluster-based subsystems are included, but experimental subsystems
are only included if the C<-usable> option is turned off.

=back

=item RETURN

Returns a reference to a hash that maps each incoming role ID to a list of
subsystem names.

    $roleHash = { $role1 => [$ss1a, $ss1b, ...],
                  $role2 => [$ss2a, $ss2b, ...],
                  ... };

=back

=cut

sub subsystems_for_role {
    # Get the parameters.
    my ($self, $args) = @_;
    # Get the sapling database.
    my $sapling = $self->{db};
    # Create the filter clause. It contains at least a role filter.
    my $filter = 'Role(id) = ?';
    # If we are NOT getting auxiliary roles, filter on the aux flag.
    if (! $args->{-aux}) {
        $filter .= " AND IsIncludedIn(auxiliary) = 0";
    }
    # There may also be subsystem filtering.
    ServerThing::AddSubsystemFilter(\$filter, $args);
    # Declare the return variable.
    my $retVal = {};
    # Get the role IDs from the parameters.
    my $ids = ServerThing::GetIdList(-ids => $args);
    foreach my $role (@$ids) {
        my @resultRows = $sapling->GetFlat("Role IsIncludedIn Subsystem",
                                          $filter, [$role], 'Subsystem(id)');
        $retVal->{$role} = \@resultRows;
    }
    # Return the result.
    return $retVal;
}

1;
