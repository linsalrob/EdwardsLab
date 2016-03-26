use strict;
use lib '/usr/localb/usr/local/lib64/perl5/auto/';
use lib '/usr/localb/usr/local/share/perl5/';
use JSON::PP;

my %tax;
open(IN, "strain_taxonomy.txt") || die "$! strain_taxonomy.txt";
while (<IN>) {
	chomp;
	# Kingdom Phylum  Class   Order   Family  Genus   Species Strain
	my @l=split /\t/;
	my $str = shift @l;
	$tax{$str}=\@l;
}
close IN;
	
my $data;

open(IN, "pairwise_ids.tsv") || die "$! pairwise_ids.tsv";
while (<IN>) {
	chomp;
	my ($str1, $str2, $pid)=split /\t/;
	next if ($str1 eq $str2);
	$str1 =~ s/^\S+\s+//;
	$str2 =~ s/^\S+\s+//;
	if ($tax{$str1} && $tax{$str2}) {
		if ($tax{$str1}->[0] eq $tax{$str2}->[0]) {push @{$data->{'kingdom'}}, $pid}
		if ($tax{$str1}->[1] eq $tax{$str2}->[1]) {push @{$data->{'phylum'}}, $pid}
		if ($tax{$str1}->[2] eq $tax{$str2}->[2]) {push @{$data->{'class'}}, $pid}
		if ($tax{$str1}->[3] eq $tax{$str2}->[3]) {push @{$data->{'order'}}, $pid}
		if ($tax{$str1}->[4] eq $tax{$str2}->[4]) {push @{$data->{'family'}}, $pid}
		if ($tax{$str1}->[5] eq $tax{$str2}->[5]) {push @{$data->{'genus'}}, $pid}
		if ($tax{$str1}->[6] eq $tax{$str2}->[6]) {push @{$data->{'species'}}, $pid}
		if ($tax{$str1}->[7] eq $tax{$str2}->[7]) {push @{$data->{'strain'}}, $pid}
	}
}
close IN;

print encode_json($data)


