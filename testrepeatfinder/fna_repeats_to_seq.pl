use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;
my $rob = new Rob;

my %opts;
getopts('f:r:v', \%opts);
unless ($opts{f} && $opts{r}) {
	die <<EOF;
	$0
	-f fasta file to read
	-r repeats file to read
	-v verbose output

	WARNING: This will probably only work for a single fasta sequence per file (not multi-fasta!!)
	WARNING: Do not use this for anything useful except a bit of debugging
EOF
}

my $fa = $rob->read_fasta($opts{f});
my @k=keys %$fa;
my $seq = $fa->{$k[0]};

open(IN, $opts{r}) || die "$! : $opts{r}";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	print "$_\t";
	$a[0] > $a[1] ? ($a[0], $a[1]) = ($a[1], $a[0]) : 1;
	print substr($seq, $a[0]-1, ($a[1] - $a[0]) - 1), "\t";
	if ($a[2] > $a[3]) {
		my $ss = substr($seq, $a[3]-1, ($a[2] - $a[3]) - 1);
		print $rob->rc($ss), "\n";
	}
	else {
		print substr($seq, $a[2]-1, ($a[3] - $a[2]) - 1), "\n";
	}

}



