use strict;
use Getopt::Std;
use Data::Dumper;


sub rc {
 my ($seq)=@_;
 $seq =~  tr/acgtrymkbdhvACGTRYMKBDHV/tgcayrkmvhdbTGCAYRKMVHDB/;
 $seq = reverse $seq;
 return $seq;
}



my %opts;
getopts('f:s:b:l:p:x:rwv', \%opts);
unless ($opts{f} && $opts{s}) {
	die <<EOF;
	$0
	-f fastq file to screen (required)
	-s SAGC ID (required)
	-b barcode file (default: barcodes.tsv)
	-r use reverse complement of the barcode
	-p percent of sequences that should have the same base at this position (default: 50)
	-l length of adapter to look for (default: 40bp)
	-x maximum number of sequences to parse (default: 1e10 ... probably the whole file!)
	-w write a long report
	-v verbose output

The barcodes.tsv file should have columns 'SAGC sample name', 'I5 Index', and 'I7 Index', but it doesn't matter what order the columns are in.

By default we just write the adapter sequence as I7_adapter\\tSequence (or I5), and this output does NOT include the barcode. If you want more output use the -w flag.

With the -w flag, we also write the barcode in lower case so you can easily distinguish it.

EOF
}

if (!$opts{p}) {$opts{p} = 50}
$opts{p} /= 100;

if (!$opts{x}) {$opts{x} = 1e10}
if (!$opts{l}) {$opts{l} = 40}

if (!$opts{b}) {$opts{b} = "barcodes.tsv"}
if (! -e $opts{b}) {die "Error barcodes file $opts{b} not found\n"}



my $sagc_column = -1;
my $i7_column = -1;
my $i5_column = -1;

my $i5; my $i7;

open(IN, $opts{b}) || die "$! : $opts{b}";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	map {$a[$_] =~ s/^\s+//; $a[$_] =~ s/\s+$//} (0..$#a);
	if ($sagc_column == -1) {
		map {
			$a[$_] eq "SAGC sample name" ? $sagc_column = $_: 1;
			$a[$_] eq "Index I7" ? $i7_column = $_: 1;
			$a[$_] eq "Index I5" ? $i5_column = $_: 1;
		} (0..$#a);
		if ($sagc_column == -1) {die "We expected to find a column in $opts{b} with the heading 'SAGC sample name in @a'\n"}
		if ($i7_column == -1) {die "We expected to find a column in $opts{b} with the heading 'I7 Index'\n"}
		if ($i5_column == -1) {die "We expected to find a column in $opts{b} with the heading 'I5 Index'\n"}
		next;
	}
	if ($a[$sagc_column] eq $opts{s}) {
		$i7 = $a[$i7_column];
		$i5 = $a[$i5_column];
	}
}
close IN;

if (!$i7 || !$i5) {die "Error: Can't find a sequence for $opts{s}. We found I7: '$i7' and I5: '$i5'\n"}

if ($opts{r}) {
	$i5 = &rc($i5);
	$i7 = &rc($i7);
}
my $seq_n = 0;
# should supplement this with fastq_stream
open(IN, "gunzip -c $opts{f} | ") || die "can't gunzip to a pipe for $opts{f}";
my @i5counts;
my @i7counts;
my $i7n = 0;
my $i5n = 0;
while (<IN>) {
	my $id = $_;
	my $seq = <IN>;
	my $discard = <IN>;
	my $qual = <IN>;
	chomp($seq);
	chomp($id);
	$seq_n++;
	last if ($seq_n > $opts{x});

	if ((my $ix = index($seq, $i7)) > -1) {
		if ($ix  < $opts{l}) {
			$opts{v} && print STDERR "Skippedi I7 match in $id because match at position $ix is less than -l option: $opts{l}\n";
		} else {
			$i7n++;
			my $start = $ix - $opts{l};
			my $end = $ix + length($i7);
			for (my $i=$start; $i < $end; $i++) {
				$i7counts[$i-$start]->{substr($seq, $i, 1)}++;
			}
		}
	}
	elsif ((my $ix = index($seq, $i5)) > -1) {
		if ($ix  < $opts{l}) {
			$opts{v} && print STDERR "Skippedi I5 match in $id because match at position $ix is less than -l option: $opts{l}\n";
		} else {
			$i5n++;
			my $start = $ix - $opts{l};
			my $end = $ix + length($i5);
			for (my $i=$start; $i < $end; $i++) {
				$i5counts[$i-$start]->{substr($seq, $i, 1)}++;
			}
		}
	}


}

if ($i7n) {
	my $adapter = "";
	for my $base (@i7counts) {
		my $tot = 0;
		map {$tot+=$base->{$_}} (qw[A C G T]);
		my $found = 0;
		if ($tot == 0) {
			$adapter .= ">";
			next;
		}
		foreach my $test (qw[A G C T]) {
			if (($base->{$test} / $tot) > $opts{p}) {
				$adapter .= $test;
				$found = 1;
				last;
			}
		}
		if (!$found) {
			$adapter .= ".";

		}
	}

	if ($adapter =~ m/([AGTC]+)$i7/) {print "$opts{f}\t$opts{s}\tI7 adapter\t$1\n"}
	if ($opts{w}) {
		print "I7 Matches. Sequence $i7. From $i7n sequences\n";
		for (my $i=0; $i < $opts{l}+length($i7); $i++) {
			if ($i % 10 > 1) {print "."}
			if ($i == 1) {print "."}
			if ($i % 10 == 0) {print $i}
		}
		print "\n";
		$adapter =~ s/$i7/lc($i7)/e;
		print "$adapter\n\n";
	}
}

if ($i5n) {
	my $adapter = "";
	for my $base (@i5counts) {
		my $tot = 0;
		map {$tot+=$base->{$_}} (qw[A C G T]);
		my $found = 0;
		if ($tot == 0) {
			$adapter .= ">";
			next;
		}
		foreach my $test (qw[A G C T]) {
			if (($base->{$test} / $tot) > $opts{p}) {
				$adapter .= $test;
				$found = 1;
				last;
			}
		}
		if (!$found) {
			$adapter .= ".";
		}
	} 
	if ($adapter =~ m/([AGTC]+)$i5/) {print "$opts{f}\t$opts{s}\tI5 adapter\t$1\n"}
	if ($opts{w}) {
		print "I5 Matches. Sequence $i5. From $i5n sequences\n";
		for (my $i=0; $i < $opts{l}+length($i5); $i++) {
			if ($i % 10 > 1) {print "."}
			if ($i == 1) {print "."}
			if ($i % 10 == 0) {print $i}
		}
		$adapter =~ s/$i5/lc($i5)/e;
		print "$adapter\n\n";
	}
}
