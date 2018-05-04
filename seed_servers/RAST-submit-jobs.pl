#!/usr/bin/env perl

use strict;
use RASTserver;
use Data::Dumper;
use Term::ReadKey;
$ENV{SAS_SERVER}="PUBSEED";
print STDERR "SAS is $ENV{SAS_SERVER}\n";


unless (defined $ARGV[0]) {
	die <<EOF 
$0 [options] <list of files to submit> 
[-v for verbose output] 
[-x Name prefix] 
[-n Organsm name (for fasta files)]
[-d domain can be V (virus) B (Bacteria) A (Archaea) E (Eukarya) only]
[-t taxonomy ID - you can either use a taxonomy id or a domain]

Optional
[-k kmerDataset to use]
[-r reannotate only (don't recall genes)]
[-u username for rast. Otherwise we will ask for it]
[-p password for rast.  Otherwise we will ask for it]


eg: $0 -v -x 'Mycophage ' *.fasta
EOF

}



my $rean=0;
my $kmerDataset;
my @files; my $prefix; my $verbose; my $orgname; my $orgdomain;
my $user; my $password; my $taxid;
while (@ARGV) {
	my $t=shift @ARGV;
	if ($t eq "-v") {$verbose = 1}
	elsif ($t eq "-x") {$prefix = shift @ARGV}
	elsif ($t eq "-n") {$orgname = shift @ARGV}
	elsif ($t eq "-u") {$user = shift @ARGV}
	elsif ($t eq "-p") {$password = shift @ARGV}
	elsif ($t eq "-d") {
		my $dom = shift;
		if (uc($dom) eq "V") {$orgdomain = "Virus"}
		elsif (uc($dom) eq "A") {$orgdomain = "Archaea"}
		elsif (uc($dom) eq "B") {$orgdomain = "Bacteria"}
		elsif (uc($dom) eq "E") {$orgdomain = "Eukarya"}
		else {die "$dom is not a valid domain choice. It must be A|B|E|V\n"}
	}
	elsif ($t eq "-t") {$taxid=shift @ARGV}
	elsif ($t eq "-r") {$rean=1}
	elsif ($t eq "-k") {$kmerDataset = shift @ARGV}
	else {push @files, $t}
}

if (!$user) {
	print "Please enter your RAST username:  ";
	$user = ReadLine(0);
	chomp $user;
}

if (!$password) {
	print "Please enter your RAST password:  ";
	ReadMode 2;
	$password = ReadLine(0);
	chomp $password;
	ReadMode 1;
	print "\n";
}


my $rast=new RASTserver($user, $password);
unless (defined $rast) {die "Can't connect ot the rast server"}



foreach my $file (@files) {
	open(IN, $file) || die "Can't open $file";
	my $data;
	if ($prefix) {$data->{"-organismName"}=$prefix}
	while (<IN>) {
		chomp;
		if (/^>/) {
			$data->{"-filetype"} = "fasta"; 
			$data->{"-organismName"} .= $file;
			$data->{"-organismName"} =~ s/.fasta//;
			$data->{"-organismName"} =~ s/.fa//;
			last;
		}
		if (/^\s*DEFINITION\s*(.*)$/) {
			$data->{"-organismName"}=$1; 
			$data->{"-filetype"}="Genbank";
		}
		if (/\/db_xref\=\"taxon:(\d+)/) {
			$data->{"-taxonomyID"}=$1;
		}
		if (/^ORIGIN/) {
			$data->{"-filetype"}="Genbank"; 
			last;
		}
	}

	# overwrite the information
	if (defined $orgname) {$data->{'-organismName'}=$orgname}
	if (defined $orgdomain) {$data->{"-domain"}=$orgdomain}
	if (defined $taxid) {$data->{"-taxonomyID"}=$taxid}

	unless (defined $data->{'-filetype'}) {die "$file does not appear to be a fasta file or a genbank file"}
	unless (defined $data->{'-organismName'}) {
		print "What is the name of the organism?:\n";
		$data->{'-organismName'}=<STDIN>;
		chomp($data->{'-organismName'});
	}



# generic stuff
	$data->{"-geneticCode"}=11;
	$data->{"-keepGeneCalls"}=0;
	$data->{"-geneCaller"}="RAST";
	$data->{"-file"}=$file;
	#$data->{"-taxonomyID"}=undef;

	if ($kmerDataset) {$data->{"-kmerDataset"} = $kmerDataset}
	if ($rean) {$data->{"-keepGeneCalls"} = 1}
	
	$verbose && (print Dumper($data)), "\n";

	my $jobH;
	eval {$jobH = $rast->submit_RAST_job($data)};
	if ($@) {
		print STDERR "Error\t$file\t$@\n";
		sleep 3;
		next;
	}
#print "Submitted the job. The returned result was: ", Dumper($jobH), "\n";
	if ($jobH->{'status'} eq "ok") {print "Submitted\t$file\t", $data->{'-organismName'}, "\t", $jobH->{'job_id'}, "\n"}
	else {print "Error\t$file\t", $jobH->{'error_message'}, "\nfor \n", Dumper($data)}
	sleep 2;
}
