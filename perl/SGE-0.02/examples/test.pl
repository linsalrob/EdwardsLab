#!/usr/bin/perl -w

use strict;
my @results;
my @fact = (1);
for (my $c=1; $c<=10; $c++) {
 for (my $i=1; $i<=6000; $i++) {
  my $a=1;
  for (my $y=1; $y<=$i; $y++) {
   $a+=factorial($y);
  }
  $results[$i]=$a;
 }
}

sub factorial {
    my $n = shift;
    return $fact[$n] if defined $fact[$n];
    $fact[$n] = $n * factorial($n - 1);
}


