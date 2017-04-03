# Schedule::SGE::Control

# POD docs

=head1 Schedule::SGE::Control

 Control jobs on the SGE queues. You should not use this method directly, rather you should use the SGE method that inherits from this, then all the methods herein are available to you.

=head1 AUTHOR

 Rob Edwards (rob@salmonella.org)
 3/24/05

=cut

package SGE::Control;
use strict;
use Exporter;

use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Schedule::SGE Exporter);
@EXPORT_OK = qw(qdel);
our $VERSION = '0.01';

=head2 qdel()

Delete all failed jobs from a queue (this must be run as the user who owns the jobs)

=cut

sub qdel {
 my ($self, $user)=@_;
 unless ($user) {$user =`whoami`; chomp($user)}
 print `qdel -u $user`;
}


1;
