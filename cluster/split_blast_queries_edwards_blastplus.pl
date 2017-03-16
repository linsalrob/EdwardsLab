#!/usr/bin/perl -w

# split a fasta file into a specified number of sub files and blast it against a database

use strict;
use POSIX;


my $usage=<<EOF;

$0 <options>

-f file to split
-n number to break into
-d destination directory (default = ".")

-p blast program (blastn, blastp, etc)
	-db blast database
-ex blast executable location (default is /usr/local/blast+/bin/)

-N job name (default is blastfile) 

-rev reverse the order of files that are submitted to the BLAST queue. (i.e. so you can run twice and start from the end backwards!)

-v verbose

Other things will be used as blast options. Unless -p && -db, will split and stop

EOF

my ($file,$dest, $no, $blastp, $blastdb, $jobname, $rev);
my ($eval, $outfmt) = (10, '6 std');
my $blastopt=" ";
my ($blastexec) = ('/usr/local/blast+/bin/');
my $verbose;


while (@ARGV) {
	my $test=shift(@ARGV);
	if ($test eq "-f") {$file=shift @ARGV}
	elsif ($test eq "-d") {$dest=shift @ARGV}
	elsif($test eq "-n") {$no=shift @ARGV}
	elsif($test eq "-p") {$blastp=shift @ARGV}
	elsif($test eq "-db") {$blastdb=shift @ARGV}
	elsif($test eq "-ex") {$blastexec=shift @ARGV; $blastexec.="/" unless (rindex($blastexec, "/") == length($blastexec)-1)}
	elsif($test eq "-N") {$jobname=shift @ARGV}
	elsif($test eq "-rev") {$rev=1}
	elsif ($test eq "-m") {print STDERR "blast+ does not use the -m option. We converted this to -outfmt 6\n"; $outfmt=6; my $i=shift @ARGV}
	elsif ($test eq "-outfmt") {$outfmt=shift @ARGV}
	elsif ($test eq "-e") {print STDERR "Adjusted -e to evalue\n"; $eval = shift @ARGV}
	elsif ($test eq "-v") {$verbose=1}
	elsif ($test eq "-evalue") {$eval = shift @ARGV}
	else {
		if ($test =~ m/\s/) {$blastopt .= " '" . $test . "' " }
		else {$blastopt .= " ". $test." "}
	}

}

die $usage unless ($file && $no);

if ($verbose) {print STDERR "BLAST OPTIONS: $blastopt\n"}

$jobname="blast$file" unless ($jobname);
if ($dest) {unless (-e $dest) {mkdir $dest, 0755}}
else {$dest="."}
$dest =~ s/\/$//;

# read the file and see how many > we have
if ($file =~ /gz$/) {open(IN, "gunzip -c $file |") || die "Can't open a pipe to $file"}
else {open(IN, $file)|| die "Can't open $file"}

my $counttags;
while (<IN>) {$counttags++ if (/^>/)}
close IN;
my $required=ceil($counttags/$no); # ceil rounds up so we should get less files than if we use int or real rounding.
#print STDERR "There are $counttags sequences in $file and we are going to write $required per file\n";


my $filecount=1;
my @sourcefiles;
if ($file =~ /gz$/) {open(IN, "gunzip -c $file |") || die "Can't open a pipe to $file"}
else {open(IN, $file)|| die "Can't open $file"}
$file =~ s/^.*\///;
open (OUT, ">$dest/$file.$filecount") || die "Can't open $dest/$file.$filecount";
push @sourcefiles, "$file.$filecount";

my $sofar;
while (my $line=<IN>) {
	if ($line =~ /^>/) {$sofar++}
	if (($line =~ /^>/) && !($sofar % $required) && (($counttags - $sofar) > 20)) {
# the last conditional is to make sure that we don't have a few sequences in a file at the end
		close OUT;
		$filecount++;
		open (OUT, ">$dest/$file.$filecount") || die "Can't open $dest/$file.$filecount";
		push @sourcefiles, "$file.$filecount";
	}
	print OUT $line;
}

@sourcefiles=reverse @sourcefiles if ($rev);


unless ($blastp && $blastdb) {die "Can't blast because either program or database were not specified"}

my $name="bl$file";
$name=substr($name, 0, 10);

my $pwd=`pwd`; chomp($pwd);

my $submitted=0;
foreach my $sf (@sourcefiles) {
	$submitted++;
	my $command = $blastexec . "$blastp -db $blastdb -query $pwd/$dest/$sf -out $pwd/$dest/$sf.$blastp -outfmt '$outfmt' -evalue $eval";
	$command .= " $blastopt";
	if ($verbose) {print STDERR "BLAST COMMAND: $command\n"}
	
	open(OUT, ">$dest/$name.$submitted.sh") || die "Can't open $dest/$name.$submitted.sh";
	print OUT "#!/bin/bash\n$command\n";
	close OUT;
	system("chmod a+x $dest/$name.$submitted.sh");



	# my $qsub= "qsub ";
	# $qsub .= " -o $dest -e $dest -cwd ";
	# $qsub .= " $pwd/$dest/$name.$submitted.sh";
	# print `$qsub`;
}

my $sgetaskidfile="$name.submitall";
my $c=0;
while (-e "$sgetaskidfile.$c.sh") {$c++}

open(STI, ">$sgetaskidfile.$c.sh") || die "can't write to $sgetaskidfile.$c.sh";
print STI "#!/bin/bash\n$dest/$name.\$SGE_TASK_ID.sh\n";
close STI;
`chmod +x $sgetaskidfile.$c.sh`;
mkdir "sge_output", 0755 unless (-e "sge_output");
my $output = join("", `qsub -S /bin/bash -cwd -e sge_output -o sge_output -t 1-$submitted:1 ./$sgetaskidfile.$c.sh`);
# output will be something like this:
# Your job-array 275368.1-92:1 ("bl6666666..submitall.0.sh") has been submitted
my $jobid="";
if ($output =~ /Your job-array (\d+)/) {
	$jobid=$1;
}




print STDERR $output, $submitted, " jobs submitted as array task\nJob ID: $jobid\n";


