use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;
my $rob = new Rob;

# make data suitable for use in sankeymatic. This just summarizes all the counts!

sub find_idx() {
	my ($headers, $value) = @_;
	for (my $i=0; $i<=$#$headers; $i++) {
		if ($headers->[$i] eq $value) {
			return $i;
		}
	}
}


my %opts;
getopts('f:m:v', \%opts);
unless ($opts{f}) {
	die <<EOF;
	$0
	-f file of count data to parse (required)
	-m mmseqs summary of Eukaryota Bacteria  Archaea  Viruses
	-v verbose output
EOF
}

my @headers;
my $counts;
my ($fqh, $fph, $shh, $nsh, $mmh);
open(IN, $opts{f}) || die "$! : $opts{f}";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	if (!@headers) {
		@headers = @a; 
		$fqh = &find_idx(\@headers, "fastq");
		$fph = &find_idx(\@headers, "fastq_fastp");
		$shh = &find_idx(\@headers, "sharks");
		$nsh = &find_idx(\@headers, "no_sharks");
		$mmh = &find_idx(\@headers, "mmseqs");
		next;
	}
	map {$counts->[$_] += @a[$_]} (1 .. $#a);
	if ($a[$fph] > $a[$fqh]) {print STDERR "$a[0]: More fastp than fastq reads\n"}
	if ($a[$shh] + $a[$nsh] > $a[$fph]) {print STDERR "$a[0]: More shark+no shark than fastp reads\n"}
	if ($a[$mmh] > $a[$nsh]) {print STDERR "$a[0]: More mmseqs reads than no sharks\n"}
}

# we start with fastq, so that should be everything
#

my $fq = $counts->[&find_idx(\@headers, "fastq")];
my $fp = $counts->[&find_idx(\@headers, "fastq_fastp")];
my $shark =  $counts->[&find_idx(\@headers, "sharks")];
my $noshark = $counts->[&find_idx(\@headers, "no_sharks")];
my $lost = $fp - ($shark + $noshark);
my $mmseq = $counts->[&find_idx(\@headers, "mmseqs")];
my $nosim = $noshark - $mmseq;


my ($euk, $bac, $arc, $vir, $oth) = (0, 0, 0, 0, 0);

if ($opts{m}) {
	open(IN, $opts{m}) || die "$! $opts{m}";
	while (<IN>) {
		chomp;
		if (m#Eukaryota (\d+) Bacteria (\d+) Archaea (\d+) Viruses (\d+) Other (\d+)#) {($euk, $bac, $arc, $vir, $oth) = ($1, $2, $3, $4, $5)}
	}
}


print "Please go to https://sankeymatic.com/build/ and paste this data\n\n";
print "fastq [", $fq - $fp, "] low quality\n";
print "fastq [", $fp, "] fastp\n";
print "fastp [", $shark, "] shark\n";
print "fastp [", $noshark, "] not shark\n";
if ($lost > 0) {
	print "fastp [", $lost, "] lost\n";
}
print "not shark [", $mmseq, "] sequence similarity\n";
print "not shark [", $nosim, "] unknown\n";

if ($opts{m}) {
	print "sequence similarity [$euk] Eukaryote\n";
	print "sequence similarity [$bac] Bacteria\n";
	print "sequence similarity [$arc] Archaea\n";
	print "sequence similarity [$vir] Virus\n";
	print "sequence similarity [$oth] Multiclass\n";
}


