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
	-d directory of mmseqs outputs [mmseqs]
	-v verbose output
EOF
}


# 0       UniRef50_L1JF06
# 1       1
# 2       0.636
# 3       0.636
# 4       0.338
# 5       131567
# 6       no rank
# 7       cellular organisms
# 8       Phytoene synthase (EC 2.5.1.32)
# 9       Metabolism  --> CLASS
# 10      Fatty Acids, Lipids, and Isoprenoids --> Level 1
# 11      Steroids and Hopanoids --> Level 2
# 12      Hopanoid biosynthesis --> Subsystem 
# 13      Phytoene synthase (EC 2.5.1.32) --> Function
# 14      0.5 --> Weighted count

# this is a REALLY memory inefficient way of doing this!
my $class; my $ac;
my $lvl1;  my $al1;
my $lvl2;  my $al2;
my $ss;    my $ass;
my $all;   my $aall;

my %allsub;
my %total;
my %sstotal;

opendir(DIR, $opts{d}) || die "$! : $opts{d}";
foreach my $sub (grep {$_ !~ /^\./} readdir(DIR)) {
	# mmseqs/SAGCFN_22_00789_S34/SAGCFN_22_00789_S34_tophit_report_subsystems.gz
	if (! -e "$opts{d}/$sub/${sub}_tophit_report_subsystems.gz") {
		print STDERR "No subsystems found for $sub\n";
		next;
	}
	open(IN, "gunzip -c $opts{d}/$sub/${sub}_tophit_report_subsystems.gz |") || die "$! opening pipe to  $opts{d}/$sub/${sub}_tophit_report_subsystems.gz";
	$allsub{$sub}=1;
	while (<IN>) {
		chomp;
		my @a=split /\t/;
		$total{$sub}+=$a[14]; ## this is the total of all reads
		next unless ($a[9]);
		$sstotal{$sub}+=$a[14]; ## this is the total of only those reads that have a subsystems match
		$class->{$sub}->{$a[9]} += $a[14];
		$ac->{$a[9]} = 1;
		$lvl1->{$sub}->{$a[10]} += $a[14];
		$al1->{$a[10]} = 1;
		$lvl2->{$sub}->{"$a[10]; $a[11]"} += $a[14];
		$al2->{"$a[10]; $a[11]"} = 1;
		$ss->{$sub}->{$a[12]} += $a[14];
		$ass->{$a[12]} = 1;
		$all->{$sub}->{"$a[9]; $a[10]; $a[11]; $a[12]"} += $a[14];
		$aall->{"$a[9]; $a[10]; $a[11]; $a[12]"} = 1;
	}
	close IN;
}


my @allsub = sort {$a cmp $b} keys %allsub;
mkdir "subsystems", 0755;

open(OUT, ">subsystems/class_raw.tsv") || die "$! subsystems/class_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ac) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $class->{$sub}->{$s} ? $class->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">subsystems/level1_raw.tsv") || die "$! subsystems/level1_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al1) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl1->{$sub}->{$s} ? $lvl1->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">subsystems/level2_raw.tsv") || die "$! subsystems/level2_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al2) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl2->{$sub}->{$s} ? $lvl2->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">subsystems/subsystems_raw.tsv") || die "$! subsystems/subsystems_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ass) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $ss->{$sub}->{$s} ? $ss->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">subsystems/all_raw.tsv") || die "$! subsystems/all_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$aall) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $all->{$sub}->{$s} ? $all->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;

## NORMALIZED COUNTS

# we divide the totals by 1e6 so that when we do the division we now are normalized per million reads
# NORMALIZED TO ALL READS THAT HAVE A MATCH

map {$total{$_} /= 1e6} keys %total;

open(OUT, ">subsystems/class_norm_all.tsv") || die "$! subsystems/class_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ac) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $class->{$sub}->{$s} ? $class->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">subsystems/level1_norm_all.tsv") || die "$! subsystems/level1_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al1) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl1->{$sub}->{$s} ? $lvl1->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">subsystems/level2_norm_all.tsv") || die "$! subsystems/level2_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al2) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl2->{$sub}->{$s} ? $lvl2->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">subsystems/subsystems_norm_all.tsv") || die "$! subsystems/subsystems_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ass) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $ss->{$sub}->{$s} ? $ss->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">subsystems/all_norm_all.tsv") || die "$! subsystems/all_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$aall) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $all->{$sub}->{$s} ? $all->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

# we divide the totals by 1e6 so that when we do the division we now are normalized per million reads
# NORMALIZED TO ALL READS WITH SUBSYSTEMS

map {$sstotal{$_} /= 1e6} keys %sstotal;

open(OUT, ">subsystems/class_norm_ss.tsv") || die "$! subsystems/class_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ac) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $class->{$sub}->{$s} ? $class->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">subsystems/level1_norm_ss.tsv") || die "$! subsystems/level1_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al1) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl1->{$sub}->{$s} ? $lvl1->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">subsystems/level2_norm_ss.tsv") || die "$! subsystems/level2_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al2) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl2->{$sub}->{$s} ? $lvl2->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">subsystems/subsystems_norm_ss.tsv") || die "$! subsystems/subsystems_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ass) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $ss->{$sub}->{$s} ? $ss->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">subsystems/all_norm_ss.tsv") || die "$! subsystems/all_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$aall) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $all->{$sub}->{$s} ? $all->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;



