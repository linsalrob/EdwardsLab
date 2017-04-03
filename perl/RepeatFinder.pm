#------------------------------------------------------------------
#
# BioPerl module Bio::Tools::RepeatFinder
#
# Cared for by Rob Edwards <redwards@utmem.edu>
#
# You may distribute this module under the same terms as perl itself
#------------------------------------------------------------------

## POD Documentation:

=head1 NAME

Bio::Tools::RepeatFinder - object to find repeats within sequences
Please note: This was never released with bioperl and now resides as RepeatFinder

=head1 SYNOPSIS

my $seqobj=Bio::Seq->new(-display_name=>"sequence", 
    -seq=>'ggagagagatgcgacgatacgagctgacgatcagcccatgcggagagagatgcgac');

my $rep=Bio::Tools::RepeatFinder->new(
 -seq=>$seqobj, -minimum=>$6, -near_to=>20);

my $newseq=$rep->joined_repeats;
$seqout->write_seq($newseq);
  

=head1 DESCRIPTION

RepeatFinder will look for direct or indirect repeats in a DNA sequence.  (Actually, it will 
find repeats in any sequence, I just never used it with proteins, so don't blame me if
something weird happens!).

There are two many options. The faster one is to return all the exact repeats. These
are repeats where each half has identical sequences. There is also an option
to join nearby repeats, returning imperfect repeats. The latter method is much
slower because it has to (a) find all the exact repeats, and then (b) compare
all of those will all the other repeats (well, almost).

I routinely use this with sequences upto about 100-200 kb, and it is bareable, but
don't try it with much longer sequences! (Or, do try it, and let me know how long
it takes).

It will return the repeats as Features and Location objects on a sequence.

Note: at the moment Indirect and Direct repeats are returned with the features on the same
strand. This is a feature because some programs (esp. artemis from Sanger) can't deal with
features on different strands being fixed. The code is there to fix this, but I use artemis
a lot, so I haven't fixed it! It could (should?) be added as a switch, I guess.

=head1 FEEDBACK

=head2 Mailing Lists 

 User feedback is an integral part of the evolution of this and other Bioperl
 modules. Send your comments and suggestions preferably to one of the Bioperl
 mailing lists. Your participation is much appreciated.

    bioperl-l@bioperl.org             - General discussion
    http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

 Report bugs to the Bioperl bug tracking system to help us keep track the bugs
 and their resolution. Bug reports can be submitted via email or the web:

     bioperl-bugs@bio.perl.org
     http://bugzilla.bioperl.org/

=head1 AUTHOR

 Rob Edwards <redwards@utmem.edu>

=head1 COPYRIGHT

 Copyright (c) 2003 Rob Edwards. All Rights Reserved.
 This module is free software; you can redistribute it and/or 
 modify it under the same terms as Perl itself.

=head1 APPENDIX

 Methods beginning with a leading underscore are considered private
 and are intended for internal use by this module. They are
 not considered part of the public interface and are described here
 for documentation purposes only.

=cut


#package Bio::Tools::RepeatFinder;
package RepeatFinder;
use Bio::Root::Root;
use Bio::Location::Split;
use Bio::Location::Simple;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Bio::Range;
use Tie::RefHash; # this is used for the _join_near_matches subroutine at the end. I can probably remove it, but I am not sure!
use strict;


use vars qw ($AUTOLOAD @RES %OK_FIELD @ISA $ID $version);

BEGIN {
 @RES=qw[verbose]; # nothing here yet, not sure what we want!

 foreach my $attr (@RES) {$OK_FIELD{$attr}++}
}


@ISA = qw(Bio::Root::Root);

#$ID = 'Bio::Tools::RepeatFinder';
$ID  = 'RepeatFinder';
$version = 1.0;


sub AUTOLOAD {
 my $self = shift;
 my $attr = $AUTOLOAD;
 $attr =~ s/.*:://;
 $self->throw("Unallowed parameter: $attr !") unless $OK_FIELD{$attr};
 $self->{$attr} = shift if @_;
 return $self->{$attr};
}


=head1 new

 Title     : new
 Purpose   : Initializes the repeat object
 Returns   : The Repeat object 
 Argument  : The DNA sequence object (required)
           : -forward=>1 for forward only, 
	   : -reverse=>1 for reverse only
	   : -newseq=>1 to return a new sequence object with the repeats, leaving the original untouched (see notes below)
	   : -left_start for the leftmost start of the region to look for repeats
	   : -left_end for the leftmost end of the region to look for repeats
	   : -right_start for the rightmost end of the region to look for repeats
	   : -right_end for the rightmost end of the region to look for repeats
	   : -minimum minimum length of repeats that are required (defaults to 7nt)
	   : -near_to distance between repeats to join as a single repeat (defaults to 10nt)
 Comments  : This is the place to start. Pass in a sequence, and you will be able to get some repeats out

=cut

sub new {
   my ($class,@args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($seq, $forward, $reverse, $leftstart, $leftend, $rightstart, $rightend, $minimum, $verbose, $near_to, $newseq) = 
     $self->_rearrange([qw(SEQ FORWARD REVERSE LEFT_START LEFT_END RIGHT_START RIGHT_END MINIMUM VERBOSE NEAR_TO NEWSEQ)],
                                                           @args);
   if (! $seq) {$self->throw("You must supply a sequence object")}
   
   
   unless ($forward || $reverse) {$forward=1; $reverse=1}
   unless ($minimum) {$minimum=7}
      
   defined $seq && $self->seq($seq);
   defined $forward && $self->set_forward($forward);
   defined $reverse && $self->set_reverse($reverse);
   defined $leftstart && $self->left_start($leftstart);
   defined $leftend && $self->left_end($leftend);
   defined $rightstart && $self->right_start($rightstart);
   defined $rightend && $self->right_end($rightend);
   defined $minimum && $self->minimum($minimum);
   defined $verbose && $self->verbose($verbose);
   defined $near_to && $self->near_to($near_to);
   defined $newseq && $self->newseq($newseq);
   
   return $self;
}

=head1 seq

 Title     : seq
 Purpose   : Get or set the sequence object
 Returns   : A Bio::Seq object
 Argument  : A Bio::Seq object
 Comments  : 

=cut

sub seq {
 my ($self, $val)=@_;
 if ($val) {$self->{seq}=$val}
 return $self->{seq};
}

=head1 set_forward

 Title     : set_forward
 Purpose   : Get or set whether forward repeats are set
 Returns   : Boolean
 Argument  : 1 to set forward repeats to be checked, -1 to unset.
 Comments  : 

=cut

sub set_forward {
 my ($self, $val)=@_;
 if ($val==1) {$self->{forward}=$val}
 elsif ($val==-1) {$self->{forward}=0}
 return $self->{forward};
}

=head1 set_reverse

 Title     : set_reverse
 Purpose   : Get or set whether reverse repeats are set
 Returns   : Boolean
 Argument  : 1 to set reverse repeats to be checked, -1 to unset.
 Comments  : 

=cut

sub set_reverse {
 my ($self, $val)=@_;
 if ($val==1) {$self->{reverse}=$val}
 elsif ($val==-1) {$self->{reverse}=0}
 return $self->{reverse};
}

=head1 newseq

 Title     : newseq
 Purpose   : Get or set whether the repeats are returned on a new sequence object
 Returns   : Boolean
 Argument  : 1 to set the sequence to a new object
 Comments  : This is not perfect, and at the moment loses a lot of sequence info, but
           : I needed it for something else I was doing at the time.

=cut

sub newseq {
 my ($self, $val)=@_;
 $self->{newseq}=$val;
 return $self->{newseq}
}

=head1 left_start

 Title     : left_start
 Purpose   : Get or set the leftmost end to begin looking for repeats
 Returns   : The position to begin looking for repeats at the left side of the sequence
 Argument  : The position
 Comments  : 

=cut

sub left_start {
 my ($self, $val)=@_;
 if ($val) {$self->{left_start}=$val}
 return $self->{left_start};
}

=head1 left_end

 Title     : left_end
 Purpose   : Get or set the end of the left side to begin looking for repeats
 Returns   : The position to end looking for repeats at the left end
 Argument  : The position
 Comments  : 

=cut

sub left_end {
 my ($self, $val)=@_;
 if ($val) {$self->{left_end}=$val}
 return $self->{left_end};
}

=head1 right_start

 Title     : right_start
 Purpose   : Get or set the rightmost end to begin looking for repeats
 Returns   : The position to begin looking for repeats at the right side of the sequence
 Argument  : The position
 Comments  : 

=cut

sub right_start {
 my ($self, $val)=@_;
 if ($val) {$self->{right_start}=$val}
 return $self->{right_start};
}

=head1 right_end

 Title     : right_end
 Purpose   : Get or set the end of the right side to begin looking for repeats
 Returns   : The position to end looking for repeats at the right end
 Argument  : The position
 Comments  : 

=cut

sub right_end {
 my ($self, $val)=@_;
 if ($val) {$self->{right_end}=$val}
 return $self->{right_end};
}

=head1 minimum

 Title     : minimum
 Purpose   : Get or set the minimum length of repeats to look for
 Returns   : The minimum length of sequence that will be checked
 Argument  : The length
 Comments  : 

=cut

sub minimum {
 my ($self, $val)=@_;
 if ($val) {$self->{minimum}=$val}
 return $self->{minimum};
}


=head1 forward_repeats

 Title     : forward_repeats
 Purpose   : Retrieve all the forward repeats
 Returns   : A hash with the key as the start location and the value as an array of the lengths
 Argument  : None
 Comments  : not implemented yet

=cut

sub forward_repeats {
 my $self=shift;
 $self->throw("forward_repeats not implemented yet");
}

=head1 reverse_repeats

 Title     : reverse_repeats
 Purpose   : Retrieve all the reverse repeats
 Returns   : A hash with the key as the start location and the value as an array of the lengths
 Argument  : None
 Comments  : not implemented yet

=cut

sub reverse_repeats {
 my $self=shift;
 $self->throw("reverse_repeats not implemented yet");
}

=head1 repeats

 Title     : repeats
 Purpose   : Retrieve a sequence object with the repeats
 Returns   : A Bio::Seq object
 Argument  : None
 Comments  : 

=cut

sub repeats {
 my $self=shift;
 if ($self->{newseq}) {
   # just add each of the features to near_to_sequence
   # for each sequence, start at the beginning and look for repeats.
   # copy the sequence so we can add repeats to it.
   $self->{repeatseq}=Bio::Seq->new(-seq=>$self->{seq}->seq);
   foreach my $feat ($self->{seq}->get_all_SeqFeatures) {$self->{repeatseq}->add_SeqFeature($feat)}
  }
  else {$self->{repeatseq}=$self->{seq}}
 $self->_find_repeats;
 return $self->{repeatseq};
}


=head1 joined_repeats

 Title     : joined_repeats
 Purpose   : Retrieve a new sequence object with the repeats where nearby repeats have been merged
 Returns   : A Bio::Seq object
 Argument  : None
 Comments  : This is a CPU intensive process because it has to look through all combinations of matches
           : and is recursive. Therefore, this is not calculated unless it is requested. It can take
	   : more than 30-40 minutes to run this process, and you might want to enable verbose to
	   : see that it is working.

=cut

sub joined_repeats {
 my $self=shift;
 # for each sequence, start at the beginning and look for repeats.
 # copy the sequence so we can add repeats to it.
 unless ($self->{fwdmatch} && $self->{revmatch}) {
  $self->{repeatseq}=Bio::Seq->new(-seq=>$self->{seq}->seq);
  foreach my $feat ($self->{seq}->get_all_SeqFeatures) {$self->{repeatseq}->add_SeqFeature($feat)}
  $self->_find_repeats;
 }
 unless ($self->{near_to_seq}) {
  if ($self->{newseq}) {
   # just add each of the features to near_to_sequence
   # for each sequence, start at the beginning and look for repeats.
   # copy the sequence so we can add repeats to it.
   $self->{near_to_seq}=Bio::Seq->new(-seq=>$self->{seq}->seq);
   foreach my $feat ($self->{seq}->get_all_SeqFeatures) {$self->{near_to_seq}->add_SeqFeature($feat)}
  }
  else {$self->{near_to_seq}=$self->{seq}}
  if ($self->{verbose}) {print STDERR "Combining forward matches within ", $self->{near_to}, "\n"}
  $self->_join_near_matches('fwd', \%{$self->{fwdmatch}});
  if ($self->{verbose}) {print STDERR "Combining reverse matches within ", $self->{near_to}, "\n"}
  $self->_join_near_matches('rev', \%{$self->{revmatch}});
 }
 return $self->{near_to_seq};
}

=head1 near_to

 Title     : near_to
 Purpose   : Get/set the near_to value
 Returns   : the near_to value
 Argument  : optional value to set
 Comments  : This sets the overlap for the repeats that are near to each other and should
           : be joined together. The default is 10 nt on EITHER side of the repeat. Note
           : it is not possible yet to specify different values for each side of the repeat,
           : and so any near_to value is applied to both sides of the repeat

=cut

sub near_to { 
 my ($self, $val)=@_;
 if ($val) {$self->{near_to}=$val}
 unless ($self->{near_to}) {$self->{near_to}=10}
 return $self->{near_to};
}


=head1 _find_repeats

 Title     : _find_repeats
 Purpose   : An internal method to figure out the repeats
 Returns   : Nothing
 Argument  : None
 Comments  : 

=cut

sub _find_repeats {
 my $self=shift;
 
 my $seq = $self->{seq}->seq; # just save us from writing the hash each time
 my $length = $self->{seq}->length;
 
 # before we carry on, need to check whether right and left ends are defined.
 # if not, set them to 1, and the length of the sequence
 
 unless ($self->{left_start}) {$self->left_start(1)}
 unless ($self->{left_end}) {$self->left_end($length)}
 unless ($self->{right_start}) {$self->right_start(1)}
 unless ($self->{right_end}) {$self->right_end($length)}
 

 my $repeatsize = $self->{minimum};
 if ($self->{verbose}) {print STDERR " $length bp\n"}


 my $repeatcount; # a counter to add to the feature
 # looking in the forward direction to begin with
 if ($self->{forward}) {
  my %fwdmatch; # forward matches
  my %match; # matches we have seen
  if ($self->verbose) {print STDERR "\nFORWARD for ", $self->{seq}->display_id, "\n"}
  my $posn=1;
  while ($posn <= ($length-$repeatsize)) {
  
   # we don't need to carry on unless this is in the left end
   unless (($posn >= $self->{left_start}) && ($posn <= $self->{left_end})) {$posn++; next}
   # just print out how far we have got
#   if ($self->verbose) {unless ($posn%(int($length/100))) {print STDERR "$posn "}}

   # now we need to get the sequence to test, and we shall start testing at 1 base futher than we start at
   my $shortseq = $self->{seq}->subseq($posn, $posn+$repeatsize);
   my $stringpos = $posn+1;

 
   # first, find where all the repeats are for this string
   my @match;
   while (($stringpos = index($seq, $shortseq, $stringpos))>-1) {
    $stringpos++; # note - bioperl uses 1 as index but index uses 0
    push @match, $stringpos;
    $stringpos++;
   }
   
   # this should be a range method, but I am getting tired and want it to work

   # is this a sub-repeat of another repeat?
   # we are going to save the repeats as
   # match{start}->{end}->{other_start}={other_end}
   foreach my $match (@match) {  
   
    #we don't need to carry on unless this is in the right end
    next unless (($match >= $self->{right_start}) && ($match <= $self->{right_end}));
    my $submatch;
    foreach my $s1 (keys %match) {
     last if ($submatch);
     foreach my $e1 (keys %{$match{$s1}}) {
      last if ($submatch);
      foreach my $s2 (keys %{$match{$s1}->{$e1}}) {
       last if ($submatch);
       my $e2=$match{$s1}->{$e1}->{$s2};
       if (($posn > $s1) && ($posn < $e1) && ($match > $s2) && ($match < $e2)) {
	$submatch=1;
       }
      }
     }
    }

    next if ($submatch);

    $match--; # note, bioperl uses 1 as index as index uses 0
    # now check for a match at posn-repeat+1
    my $longestmatch = $repeatsize+1;
    my $testlonger = index($seq,  $self->{seq}->subseq($posn, $posn+$longestmatch), $match);
    # we only want this match if it is near the original, and not others!
    unless ($testlonger == $match) {$testlonger=-1}

    # while there is still a match, increment the longest possible match
    # and test it again
    while ($testlonger > -1) {
     $longestmatch++;
     $testlonger = index($seq,  $self->{seq}->subseq($posn, $posn+$longestmatch), $match);
     unless ($testlonger == $match) {$testlonger=-1} # it is a match elsewhere.
    }
  
    $match++; # note, bioperl uses 1 as index as index uses 0
    # save the first position, the match to position and the length of the match
    # note that longest match needs to be one less than it is so that when you add
    # it to the first base position you get the last base position. I hate indices and I hate
    # that bioperl indexes on 1. That is stupid.
    $fwdmatch{$posn}->{$match}=$longestmatch-1;
     
    # save the regions that we have already seen so we don't duplicate them
    $match{$posn}->{$posn+$longestmatch}->{$match}=$match+$longestmatch-1;

    if ($self->verbose) {
     print STDERR "Added a match of ", $fwdmatch{$posn}->{$match}, " to ";
     print STDERR "position: $posn and match $match\n";
    }
    
   }
  $posn++;
  }

  foreach my $s1 (sort {$a <=> $b} keys %fwdmatch) {
   foreach my $s2 (sort {$a <=> $b} keys %{$fwdmatch{$s1}}) {
    my $location=Bio::Location::Split->new();
    $location->add_sub_Location(new Bio::Location::Simple(-start=>$s1, -end=>$s1+$fwdmatch{$s1}->{$s2}, -strand=>1));
    $location->add_sub_Location(new Bio::Location::Simple(-start=>$s2, -end=>$s2+$fwdmatch{$s1}->{$s2}, -strand=>1));
    $repeatcount++;
    my $ntmatch=$fwdmatch{$s1}->{$s2}+1;
    my $feature=Bio::SeqFeature::Generic->new(-primary => 'repeat', -source_tag   => 'repeatfinder',
      -display_name => "direct repeat $repeatcount", -tag=>{note=>"Direct repeat $repeatcount $ntmatch nt"});
    $feature->location($location);
    $self->{repeatseq}->add_SeqFeature($feature);
    #$self->{seq}->add_SeqFeature($feature);
   }
  }
 %{$self->{fwdmatch}}=%fwdmatch;
 #$self->_join_near_matches('fwd', \%fwdmatch);
 }
 

 # now looking at the reverse complements
 if ($self->{reverse}) {
  my %revmatch; # rev matches
  if ($self->verbose) {print STDERR "\nREVERSE for ", $self->{seq}->display_id, "\n"}
  my $posn=1; my %match;
  while ($posn <= ($length-$repeatsize)) {
   # we don't need to carry on unless this is in the left end
   unless (($posn >= $self->{left_start}) && ($posn <= $self->{left_end})) {$posn++; next}
   
   # just print out how far we have got
#   if ($self->verbose) {unless ($posn%(int($length/100))) {print STDERR "$posn "}}

   # now we need to get the sequence to test, and we shall start testing at 1 base futher than we start at
   my $shortseq = $self->{seq}->trunc($posn, $posn+$repeatsize)->revcom->seq;
   my $stringpos = $posn+1;
 
   # first, find where all the repeats are for this string
   my @match;
   while (($stringpos = index($seq, $shortseq, $stringpos))>-1) {
    $stringpos++; # note - bioperl uses 1 as index but index uses 0
    push @match, $stringpos;
    $stringpos++;
   }
   # is this a sub-repeat of another repeat?
   # we are going to save the repeats as
   # match{start}->{end}->{other_start}={other_end}
   foreach my $match (@match) {
   
    #we don't need to carry on unless this is in the right end
    next unless (($match >= $self->{right_start}) && ($match <= $self->{right_end}));
    
    my $submatch;
    foreach my $s1 (keys %match) {
     last if ($submatch);
     foreach my $e1 (keys %{$match{$s1}}) {
      last if ($submatch);
      foreach my $s2 (keys %{$match{$s1}->{$e1}}) {
       last if ($submatch);
       my $e2=$match{$s1}->{$e1}->{$s2};
       if (($posn >= $s1) && ($posn <= $e1) && ($match >= $s2) && ($match <= $e2)) {$submatch=1}
      }
     }
    }

    next if ($submatch);

    $match--; # note, bioperl uses 1 as index as index uses 0
    # now check for a match at posn-repeat+1
    my $longestmatch = $repeatsize+1;
    my $rmatch=$match-1; # $rmatch is because revcom matches must get smaller!
    my $testlonger = index($seq,  $self->{seq}->trunc($posn, $posn+$longestmatch)->revcom->seq, $rmatch);
    # we only want this match if it is near the original, and not others!
    
    unless ($testlonger == $rmatch) {$testlonger=-1}
    # while there is still a match, increment the longest possible match
    # and test it again
    while ($testlonger > -1) {
     $longestmatch++; $rmatch--;
     $testlonger = index($seq,  $self->{seq}->trunc($posn, $posn+$longestmatch)->revcom->seq, $rmatch);
     unless ($testlonger == $rmatch) {$testlonger=-1} # it is a match elsewhere.
    }
    
    $rmatch+=2; # note, bioperl uses 1 as index as index uses 0, and we need to counter the last --
    # save the first position, the match to position and the length of the match
    # longest match is one less because the first base. I hate indices.
    $revmatch{$posn}->{$rmatch}=$longestmatch-1;

    # save the regions that we have already seen so we don't duplicate them
    $match{$posn}->{$posn+$longestmatch}->{$rmatch}=$rmatch+$longestmatch-1;
   if ($self->verbose) {
    print STDERR "Added a match of ", $revmatch{$posn}->{$rmatch}, " to ";
    print STDERR "position: $posn and match $rmatch\n";
   }
  }
  $posn++;
  }

  foreach my $s1 (sort {$a <=> $b} keys %revmatch) {
   foreach my $s2 (sort {$a <=> $b} keys %{$revmatch{$s1}}) {
    my $location=Bio::Location::Split->new();
    $location->add_sub_Location(new Bio::Location::Simple(-start=>$s1, -end=>$s1+$revmatch{$s1}->{$s2}, -strand=>1));
    $location->add_sub_Location(new Bio::Location::Simple(-start=>$s2, -end=>$s2+$revmatch{$s1}->{$s2}, -strand=>1));
    $repeatcount++;
    my $ntmatch=$revmatch{$s1}->{$s2}+1;
    my $feature=Bio::SeqFeature::Generic->new(-primary => 'repeat', -source_tag   => 'repeatfinder',
      -display_name => "Indirect repeat $repeatcount", -tag=>{note=>"Indirect repeat $repeatcount $ntmatch nt"});
    $feature->location($location);
    $self->{repeatseq}->add_SeqFeature($feature);
    #$self->{seq}->add_SeqFeature($feature);
   }
  }
 %{$self->{revmatch}}=%revmatch;
 #$self->_join_near_matches('rev', \%revmatch);
 }
}

=head1 _join_near_matches

 Title     : _join_near_matchess
 Purpose   : An internal method to join matches that are near to each other
 Returns   : Nothing
 Argument  : None
 Comments  : "Near" is defined with by setting "near_to" in the new constructor
           : the default is 10 bp

=cut


sub _join_near_matches {
 my ($self, $direction, $match)=@_;
 unless ($self->{near_to}) {$self->{near_to}=10}
 my %match=%$match; #just dereference it so I don't forget later.....
 my %joined;
 # we are going to try and do this with Range.pm
 # Here is the theory: set up a range for each pair of starts/stops with the added near_to value
 # then find out which ranges overlap. Save those as a single location and others that don't overlap
  
 # however, we need to set up two ranges for each of the overlaps. This is not as easy as I want!
 my %range; 
 tie %range, "Tie::RefHash::Nestable";
 foreach my $s1 (sort {$a <=> $b} keys %match) {
  foreach my $s2 (sort {$a <=> $b} keys %{$match{$s1}}) {
   # start1/2 and end1/2 have the elongated starts and stops
   my ($start1, $end1)=($s1-$self->{near_to}, $s1+$match{$s1}->{$s2}+$self->{near_to});
   my ($start2, $end2)=($s2-$self->{near_to}, $s2+$match{$s1}->{$s2}+$self->{near_to});
   # now define the ranges
   my $strand;
   if ($direction eq 'fwd') {$strand=1}
   elsif ($direction eq 'rev') {$strand=-1}
   else {die "Huh, what is $direction"}
   my $range1=Bio::Range->new(-start=>$start1, -end=>$end1, -strand=>1); # the first sequence is ALWAYS on the top strand
   my $range2=Bio::Range->new(-start=>$start2, -end=>$end2, -strand=>$strand);
   $range{$range1}=$range2;
  }
 }
 my $changed=1; my $try;
 while ($changed) {
  $try++;
  my %newrange; my %changed;
  tie %newrange, "Tie::RefHash::Nestable";
  undef $changed;
  my @ranges=sort {$a->start <=> $b->start} keys %range;
  foreach my $x (0..$#ranges-1) {
   foreach my $y ($x+1..$#ranges) {
   my $rangea=$ranges[$x];
   my $rangeb=$ranges[$y];
    # having a problem where getting same start and stop all the time
    next if ($rangea->start == $rangeb->start && $rangea->end == $rangeb->end && $rangea->strand == $rangeb->strand);# the ranges are the same
     
    
    if ($rangea->overlaps($rangeb)) {
     # the ranges overlap. Do their partners?
     if ($range{$rangea}->overlaps($range{$rangeb})) {
      # yes. These two can be joined.
      #my ($starta, $enda, $stranda)=$rangea->union($rangeb); 
      #my ($startb, $endb, $strandb)=$range{$rangea}->union($range{$rangeb});
      #
      # NOTE: Docs are wrong for Bio::Range. This doesn't work as it returns a new range!
      my $range1=$rangea->union($rangeb);
      my $range2=$range{$rangea}->union($range{$rangeb});
  if ($self->verbose) {
   print STDERR "$try: RANGE A ", $rangea->start, ", ", $rangea->end, " : ", $range{$rangea}->start, ", ", $range{$rangea}->end, " ";
   print STDERR "$try: RANGE B ", $rangeb->start, ", ", $rangeb->end, " : ", $range{$rangeb}->start, ", ", $range{$rangeb}->end, " ";
   print STDERR "$try: COMBI ", $range1->start, ", ", $range1->end, " : ", $range2->start, ", ", $range2->end, "\n";
  }
      $newrange{$range1}=$range2;
      $changed{$rangea}=1;
      $changed{$rangeb}=1;
      $changed=1;
      last;
     }
    }
   }
  }
  if ($changed) {
   foreach my $range (keys %range) {
    unless (exists $changed{$range}) {$newrange{$range}=$range{$range}}
   }
   my %seen; %range=();
   foreach my $nr (keys %newrange) {
    my ($s1, $e1, $st1, $s2, $e2, $st2)=($nr->start, $nr->end, $nr->strand, $newrange{$nr}->start, $newrange{$nr}->end, $newrange{$nr}->strand);
    unless ($seen{$s1}->{$e1}->{$st1}->{$s2}->{$e2}->{$st2}) {
     $seen{$s1}->{$e1}->{$st1}->{$s2}->{$e2}->{$st2}=1; # we only want to see this once!
     $range{$nr}=$newrange{$nr}
    }
  #%range=%newrange;
   }
  }
 }
 
 my $repeatcount=1;
 my $repeattype="Direct repeat";
 if ($direction eq "rev") {$repeattype="Indirect repeat"} 
 foreach my $range (keys %range) {
  my ($s1, $e1, $str1)=($range->start, $range->end, $range->strand);
  my ($s2, $e2, $str2)=($range{$range}->start, $range{$range}->end, $range{$range}->strand);
  # now we have to correct the starts, because don't forget these ranges have -$self->{near_to}
  $s1+=$self->{near_to}; $e1-=$self->{near_to};
  $s2+=$self->{near_to}; $e2-=$self->{near_to};
  my $location=Bio::Location::Split->new();
  # note the strand=>1 should be $str1 or $str2 but most programs can't deal with repeats on different strands, so we put them on the same strand
  $location->add_sub_Location(new Bio::Location::Simple(-start=>$s1, -end=>$e1, -strand=>1));
  $location->add_sub_Location(new Bio::Location::Simple(-start=>$s2, -end=>$e2, -strand=>1)); 

  ### TO DO: 
  # add something here to check either the location or (more likely) the two sequences
  # and note whether they are perfect or imperfect repeats. it would be really cool
  # to add the number/positions of the mismatches in imperfect repeats, but this is
  # so slow anyway!
  ###
  
  $repeatcount++;
  my $ntmatch=($e1-$s1)+1;
  my $feature=Bio::SeqFeature::Generic->new(-primary => 'repeat', -source_tag   => 'repeatfinder',
    -display_name => "$repeattype $repeatcount", -tag=>{note=>"$repeattype $repeatcount $ntmatch nt"});
  $feature->location($location);
  $self->{near_to_seq}->add_SeqFeature($feature);
 }
}


=head1 to do

 Comments  : There are some things that should be added, but I haven't done so yet.
           : 1. Add whether the joined repeats are perfect or imperfect repeats
	   :    - is the sequence the same for the repeat
	   :    - how many (and where are) the mismatches
	   :    - return an alignment object?
	   : 2. Code the forward_repeats and reverse_repeats routines that are not
	   :    implemented
	   : 3. Remove the dependency on Tie::hashref
	   : 4. Make it work a lot faster

=cut
		



1;
