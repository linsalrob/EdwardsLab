=pod 

convert the output of samtools coverage to a single table.

We add an additional line at the start of the output file which is the file name, so we can easily do that
with an echo, e.g. for parsing a whole directory

=cut

use strict;


my $dir = shift || die "$0 <directory of samtools coverage output>";
opendir(DIR, $dir) || die "$! $dir";
my @allfs;
my %allrs;
my $data;
foreach my $f (grep {$_ !~ /^\./} readdir(DIR)) {
	open(IN, "$dir/$f") || die "Can't read $dir/$f";
	my $fname = <IN>;
	if ($fname =~ /^#rname/) {
		print STDERR "First line of your file is the samtools header, so using the file names\n";
		$fname = $f;
	}
	$fname =~ s#^bam/##;
	$fname =~ s/.bam$//;
	push @allfs, $fname;
	while (<IN>) {
		chomp;
		next if (/^#rname/);
		my @a=split /\t/;
		#rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
		$data->{$fname}->{$a[0]}=$a[6];
		$allrs{$a[0]}=1;
	}
	close IN;
}

my @allrs = sort {$a cmp $b} keys %allrs;
print join("\t", "", @allrs), "\n";
foreach my $f (@allfs) {
	print $f;
	map {
		$data->{$f}->{$_} ? print "\t". $data->{$f}->{$_} : print "\t0";
	} @allrs;
	print "\n";
}









