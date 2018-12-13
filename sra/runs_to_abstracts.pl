use strict;
use DBI;

my $mon = "Nov_2018";
print STDERR "Using $mon for the SQLite database\n";
my $dbfile = "/home3/redwards/SRA/SRAdb/$mon/SRAmetadb.sqlite";
unless (-e $dbfile) {
	die "Can not find SQLite file $dbfile";
}

my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","");
my $n=250;
my $inputfile = shift || die "List of SRR IDs, one per line";

my @SRR;
open(IN, $inputfile) || die "$! $inputfile";
my %aruns;
while (<IN>) {
	chomp;
	push @SRR, $_;
	if ($#SRR = $n) {
		my %runs = sql_extract(@SRR);
		map {push @{$aruns{$_}}, @{$runs{$_}}} keys %runs;
		undef @SRR;
	}
}



if (@SRR) {
	my %runs = sql_extract(@SRR);
	map {push @{$aruns{$_}}, @{$runs{$_}}} keys %runs;
}

foreach my $t (keys %aruns) {
	print "$t\t", join(",", @{$aruns{$t}}), "\n";
}

sub sql_extract() {
	my $q = join('", "', @_);
	$q = '"' . $q . '"';
	my $exc = $dbh->prepare("select r.run_accession, e.experiment_accession, s.study_accession, s.study_title, s.study_abstract from study s inner join experiment e on s.study_accession = e.study_accession inner join run r on e.experiment_accession = r.experiment_accession where r.run_accession in ($q);");
	$exc->execute || die $dbh->errstr;
	my %runs;
	while (my @r = $exc->fetchrow_array()) {
		my $run = shift @r;
		my $t=join("\t", @r);
		push @{$runs{$t}}, $run;
	}
	return %runs;
}

