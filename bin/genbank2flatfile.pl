#!/usr/bin/perl -w

# get everything out of a genbank file

use Bio::SeqIO;
use strict;
use Rob;
my $rob = new Rob;

my $usage=<<EOF;
$0 [-f print whole DNA fasta] [-v verbose] <list of genbankfiles>

EOF

die $usage unless ($ARGV[0]);


my @files = ();
my $verbose = 0;
my $print_fasta  = 0;

foreach my $t (@ARGV) {
	if ($t eq "-f") {$print_fasta = 1}
	elsif ($t eq "-v") {$verbose = 1}
	elsif (-e $t) {push @files, $t}
	else {
		print STDERR "DOn't understand command line argument $t and it is not a file. Skipped\n";
	}
}


open(TB, ">tbl") || die "can't open tbl";

my $c;
foreach my $file (@files)
{
	my $fastafile=$file;
	$fastafile =~ s/.gbk/.fasta/;
	while (-e $fastafile) {$fastafile.=".fasta"}
	my $sio=Bio::SeqIO->new(-file=>$file, -format=>'genbank');
	my $sout;
	if ($print_fasta) {
		$sout = Bio::SeqIO->new(-file=>">$fastafile", -format=>'fasta');
	}
	while (my $seq=$sio->next_seq) {
		my $seqname=$seq->display_name;
		$verbose && print STDERR "Parsing $seqname\n";
		$print_fasta && $sout->write_seq($seq);
		my @seqdata = ( $seq->display_name.".".$seq->seq_version(), $seq->length(), $seq->desc(), $seq->primary_id());
		foreach my $feature ($seq->top_SeqFeatures()) {
			
			$c++;
			my $id; # what we will call the sequence
			my ($trans, $gi, $geneid, $prod, $np);
			my $locus = "";

			eval {$locus = join " ", $feature->each_tag_value("locus_tag")};
			eval {$trans = join " ", $feature->each_tag_value("translation")};
			eval {$np = join " ", $feature->each_tag_value("protein_id")};
			if (!$np) {
				eval {$np = join " ", $feature->each_tag_value("product")};
			}
			if (!$np) {
			    $np = $locus;
			}
			if ($feature->primary_tag ne 'CDS') {
				$verbose && print STDERR "Skipped a ", $feature->primary_tag, " feature: $np\n";
				next;
			}
			if (!$np) {
				$verbose && print STDERR "No NP for $trans. Skipped\n";
				next;
			}
			elsif (!$trans) {
				$verbose && print STDERR "No translation for $np ($trans). Skipped\n";
				next;
			}
			next unless ($trans && $np);

			eval {
				foreach my $xr ($feature->each_tag_value("db_xref")) 
				{
					($xr =~ /GI/) ? ($gi = $xr) : 1;
					($xr =~ /GeneID/) ? ($geneid = $xr) : 1;
				}
			};


			eval {$prod  = join " ", $feature->each_tag_value("product")};

			my $oids = "locus:$locus"; 
			($geneid)  && ($oids.="$geneid;");
			($gi)      && ($oids.="$gi;");
			$oids =~ s/\;$//;

## columns: sequence name, fasta file for dna sequence, version, length of sequence, genome name, genome id, protein id, protein start, protein stop, protein strand, protein sequence, protein function, protein aliases

			unless ($prod)  {
				$verbose && print STDERR "No product for $np\n";
				$prod="hypothetical protein";
			}
			my $dna = $seq->subseq($feature->start, $feature->end);
			if ($feature->strand == -1) {$dna=$rob->rc($dna)}
			print TB join("\t", $seqname, $fastafile, @seqdata, $np, $feature->start, $feature->end, $feature->strand, $trans, $dna, $prod, $oids), "\n";
			# print TB join("\t", $feature->start, $feature->end, $feature->strand, $dna), "\n\n";
		}
	}
}


