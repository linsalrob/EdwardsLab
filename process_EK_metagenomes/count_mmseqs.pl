use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;
my $rob = new Rob;

my %opts;
getopts('d:v', \%opts);
unless ($opts{d}) {
	die <<EOF;
	$0
	-d directory of mmseqs output (mmseqs)
	-v verbose output
EOF
}

# mmseqs/SAGCFN_22_00789_S34/SAGCFN_22_00789_S34_tophit_aln.gz
my $euk; my $bac; my $arc; my $vir; my $total;

opendir(DIR, $opts{d}) || die "$! : $opts{d}";
foreach my $f (grep {$_ !~ /^\./} readdir(DIR)) {


	open(IN, "gunzip -c $opts{d}/$f/${f}_tophit_aln.gz |") || die "$! : gunzip -c $opts{d}/$f/${f}_tophit_aln.gz";
	my %r1; my %r2;
	while (<IN>) {
		my @a=split /\t/;
		if (substr($a[0], -1) == 1) {
			$r1{$a[0]}++;
		} elsif (substr($a[0], -1) == 2) {
			$r2{$a[0]}++;
		} else {
			print STDERR "Huh: $a[0]\n";
		}
	}
	close IN;
	
	print "$opts{d}/${f}_R1\t", scalar(keys %r1), "\n$opts{d}/${f}_R2\t", scalar(keys %r2), "\n";

	open(IN, "gunzip -c $opts{d}/$f/${f}_report.gz |") || die "$! : gunzip -c $opts{d}/$f/${f}_report.gz";
	while (<IN>) {
		chomp;
		my @a=split /\s+/;
		if ($a[$#a] eq "root") {$total += $a[1]}
		if ($a[$#a] eq "Eukaryota") {$euk += $a[1]}
		if ($a[$#a] eq "Bacteria") {$bac += $a[1]}
		if ($a[$#a] eq "Archaea") {$arc += $a[1]}
		if ($a[$#a] eq "Viruses") {$vir += $a[1]}
	}
	close IN;
	

	# 83.5653 1455816 5871    superkingdom    2759        Eukaryota
	# 2.2665  39485   7678    superkingdom    2           Bacteria
	# 0.0512  892     21      superkingdom    2157        Archaea
	# 0.0187  326     2       superkingdom    10239     Viruses

}

print STDERR "Eukaryota $euk Bacteria $bac Archaea $arc Viruses $vir Other ", ($total-($euk + $bac + $arc + $vir)), "\n";



