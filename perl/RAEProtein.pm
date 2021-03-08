# a collection of simple seed stuff for dealing with proteins and their annotations

package RAEProtein;
use strict;

=head1 new

Instantiate an RAEProtein.
my $prot = RAEProtein->new();

=cut

sub new {
	my ($class, $args)=@_;
	my $self         = {};

	return bless $self, $class;
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



=head1 strict_hypothetical

A limited subset of the is_hypothetical

Returns true is this function is part of that subset

=cut

sub strict_hypothetical {
	my ($self, $x)=@_;
	if (! $x)                             { return 1 }
	if (lc($x) eq "hypothetical protein") { return 1 }
	if (lc($x) eq "conserved protein")    { return 1 }
	if (lc($x) eq "uncharacterized protein")    { return 1 }
	if (lc($x) eq "conserved hypothetical protein") {return 1}
	if ($x =~ /hypothetical/) {return 1}
	return 0;
}




1;
