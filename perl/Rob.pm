# a collection of stuff that I want for perl

package Rob;
use strict;
use lib './';
$SIG{INT}=sub{die "$$ is done and is exiting\n"};
$SIG{CHLD}='IGNORE';

=head1 new

 instantiate

=cut

sub new {
 my ($class) = @_;
 my $self = {};
 bless $self, $class;
 return $self
}



=head1 clean_genome

 remove the crap out of genome names so you can open them with treeing programs

 my $cleangenome=clean_genome($genome);

=cut

sub clean_genome {
 my ($self, $genome)=@_;
 $genome =~ s/\'//g; $genome =~ s/\://g;$genome=~s/\(//g; $genome=~ s/\)//g;
 return $genome;
}

=head1 read_and_exec_file
  
  read and execute some code in a file. e.g. if I have a hash stored in the file, this will execute it
  We will not actually execute the code, and that way you should be able to get the answer with
  eval rob->read_exec_file($file)

=cut

sub read_exec_file {
 my ($self, $file)=@_;
 open (IN, $file) || die "Can't open $file";
 my $read;
 while (<IN>) {
  $read .= $_;
  }
  return $read;
}

=head1 argopts

My own version of Get::Opts

Pass in a REFERENCE to @ARGV, and a list of options that you want to extract. Other things on the command line are returned as
space separated.

Note that each of the objects are returned in order.

   my ($arg1, $arg2, $arg3, $rest_args)=$rob->argopts(\@ARGV, "-a", "b", "-c");


=cut

sub argopts {
 my ($self, $argv, @argu)=@_;
 return unless ($argv);
 my %need; my $i=0;
 my @found;
 foreach my $a (@argu) 
 {
  unless ($a =~ /^-/) {$a="-".$a}
  $need{$a}=$i;
  $found[$i]=undef;
  $i++;
 }
 my $rest;
 while (@$argv)
 {
  my $test=shift @$argv;
  if (defined $need{$test}) {$found[$need{$test}]=shift @$argv}
  else {$rest .= " ".$test." "}
 }
 
 return (@found, $rest);
}


=head1 backup_limited

Backup a file by incrementing the number on the end. This does not write the new file, you need to do that

e.g.

if ($rob->backup("/var/www/cgi-bin/data/$file"), 4)
{
	open(OUT, ">/var/www/cgi-bin/data/$file");
	print OUT;
	close OUT;
}
else 
{
	die "Can't back up the file!";
}

NOTE THE DIFFERENCE: This method will limit the files to 4 total. The backup() method will make an endless number of backups.

NOTE(2): The backup method is overloaded, so you can call that with an argument and it will punt to this one!

An optional mode can be provided (see the note about that in perldoc chmod and the file will have that mode

=cut

sub backup_limited()
{
	my ($self, $file, $c, $mod)=@_;
	$c--;
	eval {
		while ($c > 0)
		{
			if (-e "$file.$c")
			{
				my $n=$c+1;
				print STDERR `cp -f $file.$c $file.$n`;
				if ($mod) {chmod $mod, "$file.$n"}
			}
			$c--;
		}
		if (-e $file) {
			`cp -f $file $file.1`;
			if ($mod) {chmod $mod, "$file.1"}
		}
	};
	return 0 if ($@);
	return 1;
}




=head1 backup

Backup a file by incrementing the number on the end. This does not write the new file, you need to do that

e.g.

if ($rob->backup("/var/www/cgi-bin/data/$file"))
{
	open(OUT, ">/var/www/cgi-bin/data/$file");
	print OUT;
	close OUT;
}
else 
{
	die "Can't back up the file!";
}

NOTE: you can overload this method, and include a maximum number of copies to keep. Eg.
$rob->backup($file, 4);

=cut

sub backup()
{
	my $self=shift;
	return $self->backup_limited(@_) if ($#_ >= 1);

	my $file=shift;
        my $c=1;
	while (-e "$file.$c") {$c++}
	`cp -f $file $file.$c`;
	return 1;
}





=head1 sum
 
 Sum the values in an array

=cut

sub sum {
	my ($self, $arr)=@_;
	return 0 unless ($arr);
	return 0 unless (@$arr);
	my $sum=0;
	foreach (@$arr) {$sum+=$_};
	return $sum;
}


=head1 mean

Calculate the mean on an array
e.g. my $mean=Rob->mean(\@array);

This is the arithmetical mean (i.e. sum/n)

=cut

sub mean {
		my ($self, $arr)=@_;
 return 0 unless ($arr);
 return 0 unless (@$arr);
 my $sum;
 foreach (@$arr) {$sum+=$_};
 return $sum/(scalar @$arr);
}

=head1 median

 Calculate the median for an array
 e.g. my $median=Rob->median(\@array);
 
 This is the median, the middle number on a list. Note I am not going to return mode.

=cut

sub median {
 my ($self, $arr)=@_;
 return 0 unless ($arr);
 @$arr=sort {$a <=> $b} @$arr;
 #@$arr=sort @$arr;
 my $mid=$#$arr/2;
 if ($mid == int($mid)) {
  return $$arr[$mid];
 } else {
  return ($$arr[$mid - 0.5] + $$arr[$mid + 0.5])/2;
 }
}

=head1 min

    Return the minimum value in an array.
    e.g. my $min=Rob->min($array);

=cut

sub min {
    my ($self, $arr)=@_;
    return undef unless ($arr);
    my $min=1e99;
    map {($_ < $min) ? ($min = $_) : 1} @$arr;
    return $min;
}


=head1 max

    Return the maximum value in an array.
    e.g. my $max=Rob->max($array);

=cut

sub max {
    my ($self, $arr)=@_;
    return undef unless ($arr);
    my $max=0; 
    map {($_ > $max) ? ($max = $_) : 1} @$arr;
    return $max;
}




=head1 stdev
 
 Calculate the standard deviation on an array.

 e.g. my $stdev=Rob->stdev(\@array);
      my $stdev=Rob->stdev([1,2,3,4]);

=cut

sub stdev {        
 my ($self, $arr)=@_;
 return 0 unless ($arr);
 my ($n, $sum, $square)=(0,0,0);
 foreach (@$arr) {
  next unless ($_);
  $n++;
  $sum += $_;
  $square += $_ * $_;
 }
 
 #stdev = sqrt(((n.sum(x squared))-((sum x) squared))/n(n-1))
 unless ($n*($n-1)) {print STDERR "STDEV at div by zero as $sum, $square, $n\n"; return 0}
 my $stdev = (($n * $square)-($sum * $sum))/($n*($n-1));
 if ($stdev < 0 && $stdev > -1e-5) {return 0}
 if ($stdev < 0) {print STDERR "STdev less than zero ($stdev) as $sum, $square, $n\n"}
 $stdev=sqrt($stdev);
 return $stdev;
}  

=head1 read_fasta

Read a fasta format file and return a hash with the data, but this is not just limited to fasta sequence files - it will also handle quality scores and such.

Takes a file as argument, and an optional boolean to ensure that it is a quality file

usage: 
	my $fa=$rob->read_fasta($fastafile);
	my $fa=$rob->read_fasta("fastafile.gz"); # can also handle gzipped files
	my $qu=$rob->read_fasta($qualfile, 1);

Returns a reference to a hash - keys are IDs values are values

=cut

sub read_fasta {
 my ($self, $file, $qual)=@_;
 if ($file =~ /\.gz$/) {open(IN, "gunzip -c $file|") || die "Can't open a pipe to $file"}
 elsif ($file =~ /\.zip$/) {open(IN, "unzip -p $file|") || die "Can't open a pipe to $file"}
 else {open (IN, $file) || die "Can't open $file"}
 my %f; my $t; my $s; my $newlinewarning;
 while (<IN>) {
  if (/\r/) 
  {
  	print STDERR "The fasta file $file contains \\r new lines. It should not do this. Please complain bitterly.\n" unless ($newlinewarning);
	$newlinewarning=1;
  	s/\r/\n/g; 
	s/\n\n/\n/g;
  }
  chomp;
  if (/^>/) {
   s#^>##;
   if ($t) {
    if ($qual) {$s =~ s/\s+/ /g; $s =~ s/\s$/ /; $s =~ s/^\s+//}
    else {$s =~ s/\s+//g}
    $f{$t}=$s;
    undef $s;
   }
   $t=$_;
  }
  else {$s .= " " . $_ . " "}
 }
 if ($qual) {$s =~ s/\s+/ /g; $s =~ s/\s$/ /; $s =~ s/^\s+//}
 else {$s =~ s/\s+//g}
 $f{$t}=$s;
 close IN;
 return \%f;
}


=head1 read_fastq

Takes a fastq file and returns a hash of tuples. The first element of the tuple is the sequence, the second element is the quality score.

e.g

my $tple = $rob->read_fastq("file.fastq");
foreach my $id (keys %$tple) {
	my $seq  = $tple->{$id}->[0];
	my $qual = $tple->{$id}->[1];
}

Note that the ID has the delimiter removed. You need to add an @ at the beginning of the sequence identifier and a + at the beginning of the quality score identifier:

@QWDW5:4:85
CACGTTGAGCGAAGATATCAAGGTAGACTACTGTCG
+QWDW5:4:85
::===:999===:7667555::,2222777777666


=cut

sub read_fastq {
	my ($self, $file)=@_;
	my $l=0; my $n=-1;
	if ($file =~ /\.gz$/) {open(IN, "gunzip -c $file|") || die "Can't open a pipe to $file"}
	elsif ($file =~ /\.zip$/) {open(IN, "unzip -p $file|") || die "Can't open a pipe to $file"}
	else {open (IN, $file) || die "Can't open $file"}
	my $id; my $tple;
	while (<IN>) {
		chomp;
		$n++;
		if ($n == 0) {
			s/^@//;
			$id=$_;
		}
		elsif ($n == 1) {
			$tple->{$id}->[0]=$_;
		}
		elsif ($n == 3) {
			$tple->{$id}->[1]=$_;
			$n=-1;
		}
	}
	close IN;
	return $tple;
}

=head1 stream_fastq

Stream a fastq file. This is like read_fastq but does not send the whole file, instead it sends a stream. The default is 1GB of data, but you can set that by adding a number to the call.

usage:
while ($rob->stream_fastq("filename")) {
	## do something
}

usage to set the stream size:
while ($rob->stream_fastq($filename, 1000000)) {
	## do something
}

The return value is exactly the same as read_fastq above.


=cut

sub stream_fastq {
	my ($self, $file, $size)=@_;
	my $fh;
	unless (defined $size) {$size = 100000000}
	if ($self->{'filehandles'}->{$file}) {
		$fh = $self->{'filehandles'}->{$file};
	}
	else {
		if ($file =~ /\.gz$/) {open($fh, "gunzip -c $file|") || die "Can't open a pipe to $file"}
		elsif ($file =~ /\.zip$/) {open($fh, "unzip -p $file|") || die "Can't open a pipe to $file"}
		else {open ($fh, $file) || die "Can't open $file"}
		$self->{'filehandles'}->{$file}=$fh;
	}

	unless (defined $self->{'stream_fastq_position'}) { $self->{'stream_fastq_position'} = 0; }
	my $currpos = $self->{'stream_fastq_position'};
	my $seekpos = tell $fh;
	my $data;
	while ($fh && !eof && $seekpos < ($currpos + $size)) {
		my $id = <$fh>;
		$id =~ s/^@//;
		chomp($id);
		my $seq = <$fh>;
		chomp($seq);
		my $trash  = <$fh>;
		my $qual = <$fh>;
		chomp($qual);
		$data->{$id}=[$seq, $qual];
		$seekpos = tell $fh;
	}
	
	$self->{'stream_fastq_position'} = $seekpos;
	return $data;
}



=head1 fastq_to_qual

Conver the fastq coded quality scores to numbers

	my $qual = $rob->fastq_to_qual("::===:999===:7667555::,2222777777666");

=cut

sub fastq_to_qual {
	my ($self, $qual)=@_;
	chomp $qual;
	return "" unless ($qual);
	my $ret = "";

	foreach my $c (split //, $qual) {
		my $n = ord($c);
		$n -= 33;
		$ret .= " $n";
	}
	$ret =~ s/^ //;
	return $ret;
}

=head1 qual_to_fastq

Convert the quality score numbers to fastq codes

	my $qual = $rob->qual_to_fastq("20 30 22 45 12");

=cut

sub qual_to_fastq {
	my ($self, $qual)=@_;
	chomp $qual;
	return "" unless ($qual);

	my $ret = "";
	foreach my $c (split /\s+/, $qual) {
		$c+=33;
		$ret .= chr($c);
	}
	return $ret;
}


=head1 N50

Calculate the N50 from a fasta sequence hash. The keys should be ids and values sequences

$n50 = $rob->N50(\%hash);

=cut

sub N50 {
	my ($self, $fa)=@_;
	my %length;
	my $total;
	foreach my $k (keys %$fa) {
		my $len = length($fa->{$k});
		$length{$len}++;
		$total += $len;
	}

	my @contigsizes = sort {$b <=> $a} keys %length;
	my $currsize=0;
	while ($currsize < int($total/2)) {
		my $l = shift @contigsizes;
		$currsize += ($l * $length{$l});
	}
	return $contigsizes[0];
}


=head1 rc

Reverse complement a sequence

=cut

sub rc {
 my ($self, $seq)=@_;
 $seq =~  tr/acgtrymkbdhvACGTRYMKBDHV/tgcayrkmvhdbTGCAYRKMVHDB/;
 $seq = reverse $seq;
 return $seq;
}

=head1 rand

Randomize an array using the fisher-yates shuffle described in the perl cookbook.

=cut

sub rand {
  my ($self, $array) = @_;
  my $i;
  for ($i = @$array; --$i; ) {
   my $j = int rand ($i+1);
   next if $i == $j;
   @$array[$i,$j] = @$array[$j,$i];
  }
  return $array;
}


=head1 combined_oligos_in

Find all the oligos in a sequence within a given size range, and combine all the oligos into a single list

use:
my $arr=$rob->combined_oligos_in($seq, $min, $max);

=cut

sub combined_oligos_in {
	my ($self, $seq, $min, $max)=@_;
	my $alloligos;
	foreach my $len ($min..$max)
	{
		map {$alloligos->{$_}=1} keys %{$self->oligos_in($seq, $len)};
		my $new;
		map {$new->{$_}=1} @{$self->combine_oligos([keys %$alloligos])};
		$alloligos=$new;
	}
	return [keys %$alloligos];
}


=head1 oligos_in

Find all the oligos in a fasta sequence

Returns a hash of with the keys are all the oligos of a given size that appear once or more in any sequence. The values are the number of occurences.

Note that the oligos are returned all in uppercase.

Usage:
	my $arr=$rob->oligos_in($seq, 5);

=cut

sub oligos_in {
	my ($self, $wholeseq, $len)=@_;
	$wholeseq=uc($wholeseq);
	my %oligo;
	foreach my $seq ($wholeseq, $self->rc($wholeseq)) {
		my $posn=0;
		while($posn<=length($seq)-$len)
		{
			my $s=substr($seq, $posn, $len);
			$oligo{$s}++;
			$posn++;
		}
	}
	return \%oligo;
}


=head1 combine_oligos

Combine oligos to remove any oligos that are shorter and redundant.

For example, combine

	TTT
	TTTTA
	TTTTAGG

and just use TTTTAGG

Pass in a reference to an array of sequences, and retrieve a reference to an array of combined oligos.


=cut

sub combine_oligos {
	my ($self, $arr)=@_;

	my @oligos=sort {length($b) <=> length($a)} @$arr;

	my %all;
	foreach my $o (@oligos)
	{
		my $seen=-1;
		foreach my $k (keys %all)
		{
			$seen=index($k, $o);
			last if ($seen > -1);
		}
		if ($seen == -1) {$all{$o}=1}
	}
	return [keys %all];
}

=head1 commify

Put commas in numbers. I think this comes straight from the perl cookbook and is very useful for nice displays

=cut

sub commify {
    my($self,$n) = @_;
    my(@n) = ();
    my($i);

    for ($i = (length($n) - 3); ($i > 0); $i -= 3)
    {
        unshift(@n,",",substr($n,$i,3));
    }
    unshift(@n,substr($n,0,$i+3));
    return join("",@n);
}


=head1 Tm

Please Note: This is stolen 100% from BioPerl, but since I wrote that bioperl module it is not plaigarism :)
You should preferentially use the bp module for this, but I don't want to since one system I am working on
doesn't have it at the moment and this is the quickest way to get it working.

Title   : Tm()
	Usage   : $tm = $primer->Tm(-salt=>'0.05', -oligo=>'0.0000001', -seq=>'ATTACTTATATCA')
	Function: Calculates and returns the Tm (melting temperature) of the primer
	Returns : A scalar containing the Tm.
	Args    : -salt set the Na+ concentration on which to base the calculation (default=0.05 molar).
	: -oligo set the oligo concentration on which to base the calculation (default=0.00000025 molar).
	: -seq the sequence. This is the only required argument
	Notes   : Calculation of Tm as per Allawi et. al Biochemistry 1997 36:10581-10594.  Also see
	documentation at http://biotools.idtdna.com/analyzer/ as they use this formula and
	have a couple nice help pages.  These Tm values will be about are about 0.5-3 degrees
	off from those of the idtdna web tool.  I don't know why.

	This was suggested by Barry Moore (thanks!). See the discussion on the bioperl-l
	with the subject "Bio::SeqFeature::Primer Calculating the PrimerTM"

=cut

sub Tm  {
	my ($self, %args) = @_;
	my $salt_conc = 0.05; #salt concentration (molar units)
		my $oligo_conc = 0.00000025; #oligo concentration (molar units)
		if ($args{'-salt'}) {$salt_conc = $args{'-salt'}} #accept object defined salt concentration
			if ($args{'-oligo'}) {$oligo_conc = $args{'-oligo'}} #accept object defined oligo concentration
				my $sequence;
	if ($args{'-seq'}) {$sequence=uc $args{'-seq'}}
	unless ($sequence) {print STDERR "No sequence provided for Tm calculation. returned 0\n";  return 0}
	my $length = length($sequence);
	my @dinucleotides;
	my $enthalpy;
	my $entropy;
#Break sequence string into an array of all possible dinucleotides
	while ($sequence =~ /(.)(?=(.))/g) {
		push @dinucleotides, $1.$2;
	}
#Build a hash with the thermodynamic values
	my %thermo_values = ('AA' => {'enthalpy' => -7.9,
			'entropy'  => -22.2},
			'AC' => {'enthalpy' => -8.4,
			'entropy'  => -22.4},
			'AG' => {'enthalpy' => -7.8,
			'entropy'  => -21},
			'AT' => {'enthalpy' => -7.2,
			'entropy'  => -20.4},
			'CA' => {'enthalpy' => -8.5,
			'entropy'  => -22.7},
			'CC' => {'enthalpy' => -8,
			'entropy'  => -19.9},
			'CG' => {'enthalpy' => -10.6,
			'entropy'  => -27.2},
			'CT' => {'enthalpy' => -7.8,
			'entropy'  => -21},
			'GA' => {'enthalpy' => -8.2,
			'entropy'  => -22.2},
			'GC' => {'enthalpy' => -9.8,
			'entropy'  => -24.4},
			'GG' => {'enthalpy' => -8,
			'entropy'  => -19.9},
			'GT' => {'enthalpy' => -8.4,
				'entropy'  => -22.4},
			'TA' => {'enthalpy' => -7.2,
				'entropy'  => -21.3},
			'TC' => {'enthalpy' => -8.2,
				'entropy'  => -22.2},
			'TG' => {'enthalpy' => -8.5,
				'entropy'  => -22.7},
			'TT' => {'enthalpy' => -7.9,
				'entropy'  => -22.2},
			'A' =>  {'enthalpy' => 2.3,
				'entropy'  => 4.1},
			'C' =>  {'enthalpy' => 0.1,
				'entropy'  => -2.8},
			'G' =>  {'enthalpy' => 0.1,
				'entropy'  => -2.8},
			'T' =>  {'enthalpy' => 2.3,
				'entropy'  => 4.1}
	);
#Loop through dinucleotides and calculate cumulative enthalpy and entropy values
	for (@dinucleotides) {
		unless (defined $thermo_values{$_}{enthalpy}) {print STDERR "No enthalpy for |$_|\n"; $thermo_values{$_}{enthalpy}=0}
		unless (defined $thermo_values{$_}{entropy}) {print STDERR "No entropy for |$_|\n"; $thermo_values{$_}{entropy}=0}
		$enthalpy += $thermo_values{$_}{enthalpy};
		$entropy += $thermo_values{$_}{entropy};
	}
#Account for initiation parameters
	$enthalpy += $thermo_values{substr($sequence, 0, 1)}{enthalpy};
	$entropy += $thermo_values{substr($sequence, 0, 1)}{entropy};
	$enthalpy += $thermo_values{substr($sequence, -1, 1)}{enthalpy};
	$entropy += $thermo_values{substr($sequence, -1, 1)}{entropy};
#Symmetry correction
	$entropy -= 1.4;
	my $r = 1.987; #molar gas constant
		my $tm = ($enthalpy * 1000 / ($entropy + ($r * log($oligo_conc))) - 273.15 + (12* (log($salt_conc)/log(10))));
	return $tm;
}



=head2

Remove non standard ascii characters from a string. This just strips out anything that is not easy to see in a terminal (thus shouldn't exist), and replaces them with a space.

use: $string = $rae->ascii_clean($string);

=cut

sub ascii_clean {
    my ($self, $str) = @_;
    return unless ($str);
    my %goodChrs = map {$_=>1} (9,10,13,32..127);
    $str =~ s/(.)/$goodChrs{ord($1)} ? $1 : ' '/eg;
    return $str;
}


=head2 iupac_code

Define the genetic code as a hash of iupac codes. Much of this is work by Gary Olsen, University of Illinois, Urbana Champaign. Thanks, Gary.

To use, sort your sequence and then you can get the iupac code

=cut

sub iupac_code {
	my $self = shift;
%{$self->{DNA_letter_can_be}} = (
		A => ["A"],                 a => ["a"],
		B => ["C", "G", "T"],       b => ["c", "g", "t"],
		C => ["C"],                 c => ["c"],
		D => ["A", "G", "T"],       d => ["a", "g", "t"],
		G => ["G"],                 g => ["g"],
		H => ["A", "C", "T"],       h => ["a", "c", "t"],
		K => ["G", "T"],            k => ["g", "t"],
		M => ["A", "C"],            m => ["a", "c"],
		N => ["A", "C", "G", "T"],  n => ["a", "c", "g", "t"],
		R => ["A", "G"],            r => ["a", "g"],
		S => ["C", "G"],            s => ["c", "g"],
		T => ["T"],                 t => ["t"],
		V => ["A", "C", "G"],       v => ["a", "c", "g"],
		W => ["A", "T"],            w => ["a", "t"],
		Y => ["C", "T"],            y => ["c", "t"]
		);
map {$self->{iupac}->{join("", @{$self->{DNA_letter_can_be}->{$_}})}=$_} keys %{$self->{DNA_letter_can_be}};

return $self->{iupac};
}

=head2 iupac

Convert an array of sequences to their relevant IUPAC codes. Sequences can be in any order/case but code will be upper case

my $iupac = $rob->iupac(\@seq);

=cut

sub iupac {
	my ($self, $array) = @_ ;
	my %array;
	unless (defined $self->{iupac}) {$self->iupac_code()}

	foreach my $base (@$array) {
		next if (!$base || $base eq " " || $base eq "-");
		$base=uc($base);
		if ($base eq "A" || $base eq "T" || $base eq "G" || $base eq "C") {$array{$base}=1}
		elsif ($self->{DNA_letter_can_be}->{$base}) {
			map {$array{$_}=1} @{$self->{DNA_letter_can_be}->{$base}};
		}
		else {
			die "Unknown base in alignment file: $base\n";
		}
	}

	my $code= join("", sort {$a cmp $b} grep {m/[GATC]+/} keys %array);
	return $self->{iupac}->{$code} if (defined $self->{iupac}->{$code});
	die "No IUPAC code for $code in ", Dumper(\%array);
}

=head2 Shannon's complexity calculation

Calculate the complexity of a given piece of DNA for a certain oligo length. Assumes DNA sequences!

my $hprime = $rob->shannon($dna, $wordlength); 

$wordlength must be >= 1

$hprime is in nats

=cut

sub shannon {
	my ($self, $dna, $word)=@_;
	die 'Must specify a word length >= 1' if (!defined $word || $word < 1);
	my %count;
	my $posn=0;
	while ($posn < length($dna)-$word) {
		$count{substr($dna, $posn, $word)}++;
		$posn++;
	}
	my $hp=0;
	my $poss = length($dna);
	map {$hp += ($count{$_}/$poss) * log($count{$_}/$poss)} keys %count;
	return -$hp;
}


=head1 background_job 

 Fork a child and run the job in the background. Do not wait for the child to finish.
 This is used on the submit pages to run jobs in the background.

 If the working directory is supplied we will chdir there first

 eg my $pid=$jc->background_job($command, $wd);

=cut

sub background_job {
	my ($self, $cmd, $wd)=@_;
	return unless ($cmd);
	chdir($wd) if ($wd && -d $wd);
	my $child;
	unless ($child=fork) {
		# this is what the child does
		die "cannot fork because $!" unless defined $child;
		# redirect STDOUT/STDERR/STDIN
		close(STDOUT);
		close(STDERR);
		open STDIN, "</dev/null";

		exec "$cmd > out.txt 2>out.err";
		exit;
	}

	return $child;
}
 
  

1;

