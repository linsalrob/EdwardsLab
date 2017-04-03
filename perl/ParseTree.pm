#!/usr/bin/perl

package ParseTree;
use strict;
=pod

=head1

A module written by Rob to parse a tree and return a hash. The hash has two entries, "left", and "right". This is based on code from perlmonks

http://www.perlmonks.org/?node_id=717769

=cut

=head1 new

 instantiate

=cut


sub new {
	my ($class)=@_;
	my $self={};
	$self->{'n'}=0;
	return bless $self, $class;
}



=head1 parse

Parse a tree. Return a hash

my $hash = ParseTree->parse($tree)
NOTE:
I added a method that preserves the distance, but changes the node to an array of three-ples: ["node", "name", "distance"]


=cut

sub parse {
	my ($self, $tree, $dist)=@_;
	
	# remove the trailing ;
	$tree =~ s/\;\s*$//;

	if ($dist) {
		return $self->FromStringWithDistance($tree);
	} else {
		return $self->FromString($tree);
	}

}

=head1 as_string

Returns a string object suitable for printing out, etc

my $string = $parse->as_string($tree);

=cut

sub as_string {
	my ($self, $tree)=@_;
	$self->{'as_string'}="";
	$self->write_tree($tree);
	return $self->{'as_string'};
}
	



sub FromString {
	my ($self, $in) = @_;
	# print STDERR "From String is deprecated. Please use From String with distance, but this will change the nodes to be tuples of node, distance\n";
# check to see if it has parens and a comma
	if (index($in, ",") == -1) {
		# if (rindex($in, ")") > 0 && rindex($in, ")") < rindex($in, ":")) {
		#	return (substr($in, rindex($in, ":"), length($in)), substr($in, 0, rindex($in, ")")-1));
		#}
		
		if (index($in, ":") != -1) {
			return substr($in, 0, index($in, ":"));
		}


		return $in;
	}

	my $p=0;
	if (index($in, "(") != -1) {
		my $depth = 0;
		for (my $i=0; $i <  length($in); $i++) { ## find split point
			my $c = substr $in, $i, 1;
			++$depth if $c eq '(';
			--$depth if $c eq ')';
			$p = $i and last if $c eq ',' and $depth == 1;
		}
	}
	else {
		$p=index($in, ",");
	}

	my $start = (index($in, "(") == 0) ? 1 : 0;
	#my $left = substr $in, 1, $p-1;
	#my $right = substr $in, $p+1, length($in);
	#print STDERR "Split ", $self->{n}++, "\t$left :: $right\n";
	return [$self->FromString( substr $in, $start, $p-1 ), $self->FromString( substr $in, $p+1, length($in) )];
}

sub FromStringWithDistance {
	my ($self, $in)=@_;

	# do we have a distance. If so we want to trim that off and save it for later
	my $distance="";
	my $rparen=rindex($in, ")");
	my $rcolon=rindex($in, ":");
	if ($rparen < $rcolon) {
		$distance = substr($in, $rcolon+1);
		$in = substr($in, 0, $rcolon); # adjust the string to exclude the distance
	}

	# if we don't have a comma we are just the node
	if (index($in, ",") == -1) {
		return ["node", $in, $distance];
	}
	
	my $p=0;
	# find the split point where we are at one level above nothing!
	if (index($in, "(") != -1) {
		my $depth = 0; # the depth we are in the tree
		for (my $i=0; $i <  length($in); $i++) {
			my $c = substr($in, $i, 1);
			++$depth if ($c eq '(');
			--$depth if ($c eq ')');
			if ($c eq ',' && $depth == 1) {
				$p = $i;
				#print STDERR "Set ptr to $p\n";
				last;
			}
		}
	}
	else {
		$p=index($in, ",");
	}

	# now we have the split points we need to figure out if the left and right halves begin and end with parens, and ignore them if they do
	my ($start, $endlength)=(0, length($in)-$p-1);
	if (index($in, "(") == 0 && rindex($in, ")") == length($in)-1) {
		$start=1; $endlength--;
	}
	
	# print STDERR substr($in, $start, $p-1 ), " *** ", substr($in, $p+1, $endlength), "\n\n";
	return [$self->FromStringWithDistance(substr($in, $start, $p-1 )), $self->FromStringWithDistance(substr($in, $p+1, $endlength)), $distance];
}





sub write_tree {
	my ($self, $tree)=@_;
		
	if (ref($tree) eq "ARRAY") {
		if (!defined $tree->[0]) {
			$self->{'as_string'} .= ":" . $tree->[2];
			if (ref($tree->[1]) eq "ARRAY") {$self->write_tree($tree->[1])}
			else {$self->{'as_string'} .= ":" . $tree->[1]}
			return;
		}
		if ($tree->[0] eq "node") {
			my $n=shift @$tree;
			print STDERR "JOINING ", join(";  ", @$tree), "\n";
			$self->{'as_string'} .= join(":", @$tree);
			return;
		}
		$self->{'as_string'} .= "(";
		$self->write_tree($tree->[0]);
		$self->{'as_string'} .= ",";
		$self->write_tree($tree->[1]);
		$self->{'as_string'} .= ")";
	}
	else {
		$self->{'as_string'} .= $tree;
	}
}





sub count_nodes {
	my ($self, $node)=@_;
	if (ref($node) eq "ARRAY") {
		my $cn;
		($node, $cn) = $self->count_children($node);
		push @$node, $cn;
		$node->[0] = $self->count_nodes($node->[0]);
		$node->[1] = $self->count_nodes($node->[1]);
		return $node;
	}
	else {
		return $node;
	}
}

=head1 number_of_children

Returns the number of children at any given point in the tree

=cut

sub number_of_children {
	my ($self, $tree)=@_;
	my ($t, $n)=$self->count_children($tree);
	return $n;
}

sub count_children {
	my ($self, $ref, $count)=@_;
	if (ref($ref) eq "ARRAY") {
		($ref->[0], $count) = $self->count_children($ref->[0], $count);
		($ref->[1], $count) = $self->count_children($ref->[1], $count);
		return ($ref, $count);
	}
	else {
		$count++;
		return ($ref, $count);
	}
}


=head1 rename_tree_leaves

Takes a tree object and a hash of "current name" => "new name" and renames all leaves

my $newtreeObj = $parse->rename_tree_leaves($tree, $names);

You probably want this one:

my $string = $parse->as_string($parse->rename_tree_leaves($tree, $names));
print "$string\n";

=cut

sub rename_tree_leaves {
	my ($self, $tree, $names) = @_;
	if (ref($tree) eq "ARRAY") {
		if ($tree->[0] eq "node") {
			if ($names->{$tree->[1]}) {
				$tree->[1] = $names->{$tree->[1]};
			}
			return $tree;
		}
		$tree->[0] = $self->rename_tree_leaves($tree->[0], $names);
		$tree->[1] = $self->rename_tree_leaves($tree->[1], $names);
		return $tree;
	}
	else {
		if ($names->{$tree}) {$tree=$names->{$tree}}
		return $tree;
	}
}



1;
