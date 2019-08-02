#!/usr/bin/perl -w

# authors: Adam Brenner <aebrenne@uci.edu>, http://www.ics.uci.edu/~aebrenne/
# and      Harry Managalam <harry.mangalam@uci.edu>

use strict;
use vars qw($OUTFILE $MAXSIZE $MAXAGE $BLOB $range $last_write_time $current_time 
	$age $waitcnt $stillclosed $random_amt         
);
# usage: someapp --opt1  --opt2 | zotkill.pl /path/to/shared/filename

$OUTFILE = shift; # this ingests the '/path/to/shared/filename'
$MAXSIZE = 1024;  # size of the output buffer
$MAXAGE = 2; # max age in seconds that the buffer is allowed to be unwritten
$BLOB = ""; # the buffer
$range = 100;  # 0-99  # the seed for the rng
$last_write_time = time(); # init the var to the # od sec since the epoch

if (!defined $OUTFILE || $OUTFILE eq "") {
	print STDERR << "ENDHELP";
  
  Usage: someapp --opt1  --opt2 | zotkill.pl /path/to/shared/filename
	
  zotkill.pl is used to trap the outputs of multiple simultaneous processes 
  and feed them into the same file, reducing the number of files that it takes 
  to accumulate multiple outputs.  We call this the Zillions of Tiny (ZOT) 
  file problem, hence the name of this utility.  The problem and its 
  solution are described here: <http://moo.nac.uci.edu/~hjm/Job.Array.ZOT.html
  This is zotkill version 1.0
  
ENDHELP
exit(0);
}

while (<>) { # while there's still stuff coming in on STDIN
    $BLOB .= $_; # just add STDIN to $BLOB
     $current_time = time(); # take the time
    $age = $current_time - $last_write_time;  # calc the age since last write
    if (length($BLOB) > $MAXSIZE || $age > $MAXAGE) {
    	file_still_locked($OUTFILE);  # acquire file lock
        print OUT "$BLOB";  # using flocks.
        print STDERR "w"; # indicates successful writes.
        close OUT;
        $BLOB = ""; # reset 
        $last_write_time = time(); #reset
    } 
}

file_still_locked($OUTFILE);
print OUT "$BLOB"; # dump the last bits
print STDERR "+"; # note it in the STDERR
close OUT;  # close and unlock

sub file_still_locked {
   $OUTFILE = shift; # grab the filename
   $waitcnt = 0;
   $stillclosed = 1;
   while ($stillclosed){
      $waitcnt++;
      $stillclosed = 0;
      $random_amt = rand($range)/100;
      sleep($random_amt);
      open OUT, ">> $OUTFILE" or $stillclosed = 1;
      if (!$stillclosed){ flock(OUT, 2); }
      else {print STDERR "w ";} # string of 'wwwwww ' to indicate waits
      if ($waitcnt > 10000) { print "Dying!!\n"; die "\n\nI've been waiting too long! Dying.....\n\n"}
    }
}

