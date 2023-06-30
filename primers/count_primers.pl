use strict;
use Getopt::Std;
use Data::Dumper;
use Rob;
my $rob = new Rob;

# given a list of directories, read all the files and count the primer occurrences

my $usage ="$0 <I7 primer> <I5 primer> fastq/ fastq/ fastq/";

unless ($#ARGV > 1) {die $usage}

my $i7 = shift || die $usage;
my $i5 = shift || die $usage;
my $i5r = $rob->rc($i5);
my $i7r = $rob->rc($i7);

print "Looking for I7: $i7 and I5: $i5 in @ARGV\n";
print "File\tI7\tI5\tI7rc\tI5rc\tn\n";
my $percs;

foreach my $dir (@ARGV) {
	opendir(DIR, $dir) || die "can't open $dir";
	foreach my $fq (grep {/fastq/} readdir(DIR)) {
		my ($i5n, $i7n, $i5rn, $i7rn, $n) = (0, 0, 0, 0, 0);
		foreach my $seq ($rob->stream_fastq("$dir/$fq")) {
			foreach my $id (keys %$seq) {
				if (index($seq->{$id}->[0], $i5) > -1) {$i5n++}
				if (index($seq->{$id}->[0], $i7) > -1) {$i7n++}
				if (index($seq->{$id}->[0], $i5r) > -1) {$i5rn++}
				if (index($seq->{$id}->[0], $i7r) > -1) {$i7rn++}
				$n++;
			}
		}
		print "$fq\t$i7n\t$i5n\t$i7rn\t$i5rn\t$n\n";
		push @{$percs->[0]}, $i5n/$n * 100;
		push @{$percs->[1]}, $i7n/$n * 100;
		push @{$percs->[2]}, $i5rn/$n * 100;
		push @{$percs->[3]}, $i7rn/$n * 100;

	}
}

$rob->sum($percs->[1]) && print "I7 ", $rob->mean($percs->[1]), "% +/- ", $rob->stdev($percs->[1]), "\n";
$rob->sum($percs->[0]) && print "I5 ", $rob->mean($percs->[0]), "% +/- ", $rob->stdev($percs->[0]), "\n";
$rob->sum($percs->[3]) && print "I7 rc ", $rob->mean($percs->[3]), "% +/- ", $rob->stdev($percs->[3]), "\n";
$rob->sum($percs->[2]) && print "I5 rc ", $rob->mean($percs->[2]), "% +/- ", $rob->stdev($percs->[2]), "\n";


