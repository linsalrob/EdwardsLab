#!perl

# This is a package that I wrote with Junior as some contract work. This is not released in cpan, but probably could be at some point.
# it is basically a wrapper around GD.pm that provides tracks and similar things. It is pretty darn good!

# If someone wanted a fun project, they could take this, edit it and add help documents, rewrite parts of it and put it in cpan.

use vars qw/$VERSION/;
BEGIN {
  $VERSION = do { my @r = (q$Revision: 1.1.1.1 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
package OGD::Object;
use strict;
use vars '$AUTOLOAD';

sub translate_coords {
  my ($self, @ac) = @_;
  my @own_coord = $self->coords;
  my @c = @own_coord;
  foreach my $i (0 .. $#ac) {
    my $frame_size = $i % 2 ? ($own_coord[3] - $own_coord[1]) : ($own_coord[2] - $own_coord[0]);
    my $frame_origin = $i % 2 ? $own_coord[1] : $own_coord[0];
    if ( $ac[$i] =~ /^(\.)?([+-]?\d+)p$/i) {
      # value expressed in pixels
      my ( $relative, $val ) = ($1, $2);
      if ( $relative ) {
	# It's relative to the previous point.
	# In case of rectangles, it's relative to the opposite corner
	$c[$i] = $c[$i-2] + $val;
      }
      else {
	$c[$i] = $val;
      }
    }
    else {
      # this coordinate is fractional
      my ( $relative, $val );
      if ( $ac[$i] =~ /^(\.)?([+-]?\d+)$/) {
	# it is also relative
	( $relative, $val ) = ($1, $2);
      }
      else {
	$val = $ac[$i];
      }
	
      if ( $relative ) {
	$c[$i] = $c[$i-2] + $val*$frame_size;
      }
      else {
	$c[$i] = $frame_origin + $val*$frame_size;
      }
    }
    $c[$i] = int $c[$i];
  }

  return @c;
}

sub line {
  my ($self, %arg) = @_;
  my @coord = $self->translate_coords(@{$arg{-coords}});
  my $color = $arg{-fill} || 'black';
  my $width = $arg{-width} || 1;
  my $arrow = $arg{-arrow} || 'none';
  my $arrowshape = $arg{-arrowshape} || [8, 10, 3];

  my @c = $self->coords;
  my @point;   # An array containing pairs of x- and y-coords of all
               # points in line, stored as anonymous arrays. If the
               # line has arrowheads then the first and last points #
               # are adjusted to refer to the necks of the arrowheads
               # rather than # their tips.

  foreach my $i ( 0 .. @coord/2 - 1 ) {
    $point[$i] = [@coord[2*$i, 2*$i+1]];
  }

  if ( $arrow ne 'none' ) {
    my @poly;  # the polygon for arrowhead
    my $fracHeight;                  # Line width as fraction of
                                     # arrowhead width.

    my $backup;                      # Distance to backup end points
                                     # so the line ends in the middle
                                     # of the arrowhead.

    my ( $vertX, $vertY );           # Position of arrowhead vertex
    my ( $shapeA, $shapeB, $shapeC); # Adjusted coordinates (see
                                     # explanation below)


    # The code below makes a tiny increase in the shape parameters for
    # the line.  This is a bit of a hack, but it seems to result in
    # displays that more closely approximate the specified parameters.
    # Without the adjustment, the arrows come out smaller than
    # expected.
    $shapeA = $arrowshape->[0] + 0.001;
    $shapeB = $arrowshape->[1] + 0.001;
    $shapeC = $arrowshape->[2] + $width/2.0 + 0.001;

    # If there's an arrowhead on the first point of the line, compute
    # its polygon and adjust the first point of the line so that the
    # line doesn't stick out past the leading edge of the arrowhead
    $fracHeight = ($width/2.0)/$shapeC;
    $backup = $fracHeight*$shapeB + $shapeA*(1.0 - $fracHeight)/2.0;

    if ( $arrow eq 'first' || $arrow eq 'both' ) {
      $poly[0] = $poly[10] = $point[0]->[0];
      $poly[1] = $poly[11] = $point[0]->[1];

      my $dx = $poly[0] - $point[1]->[0];
      my $dy = $poly[1] - $point[1]->[1];
      my $length = sqrt($dx*$dx + $dy*$dy);

      my ($sinTheta, $cosTheta);
      if ($length == 0) {
        $sinTheta = $cosTheta = 0.0;
      } else {
        $sinTheta = $dy/$length;
        $cosTheta = $dx/$length;
      }

      $vertX = $poly[0] - $shapeA*$cosTheta;
      $vertY = $poly[1] - $shapeA*$sinTheta;
      my $temp = $shapeC*$sinTheta;
      $poly[2] = $poly[0] - $shapeB*$cosTheta + $temp;
      $poly[8] = $poly[2] - 2*$temp;
      $temp = $shapeC*$cosTheta;
      $poly[3] = $poly[1] - $shapeB*$sinTheta - $temp;
      $poly[9] = $poly[3] + 2*$temp;
      $poly[4] = $poly[2]*$fracHeight + $vertX*(1.0-$fracHeight);
      $poly[5] = $poly[3]*$fracHeight + $vertY*(1.0-$fracHeight);
      $poly[6] = $poly[8]*$fracHeight + $vertX*(1.0-$fracHeight);
      $poly[7] = $poly[9]*$fracHeight + $vertY*(1.0-$fracHeight);

      # Polygon done.  Now move the first point towards the second so
      # that the corners at the end of the line are inside the
      # arrowhead
      $point[0]->[0] = $poly[0] - $backup*$cosTheta;
      $point[0]->[1] = $poly[1] - $backup*$sinTheta;

      # round up the coords and indicate that the values are in pixels
      @poly = map { int($_) . 'p' } @poly;

      # draw the arrow
      $self->gd_polygon('filled', -coords => \@poly, -color => $color);
    }

    # Similar arrowhead calculation for the last point of the line
    if ( $arrow eq 'last' || $arrow eq 'both' ) {
      $poly[0] = $poly[10] = $point[-1]->[0];
      $poly[1] = $poly[11] = $point[-1]->[1];
      my $dx = $poly[0] - $point[-2]->[0];
      my $dy = $poly[1] - $point[-2]->[1];
      my $length = sqrt($dx*$dx + $dy*$dy);

      my ($sinTheta, $cosTheta);
      if ($length == 0) {
	$sinTheta = $cosTheta = 0.0;
      } else {
	$sinTheta = $dy/$length;
	$cosTheta = $dx/$length;
      }
      $vertX = $poly[0] - $shapeA*$cosTheta;
      $vertY = $poly[1] - $shapeA*$sinTheta;
      my $temp = $shapeC*$sinTheta;
      $poly[2] = $poly[0] - $shapeB*$cosTheta + $temp;
      $poly[8] = $poly[2] - 2*$temp;
      $temp = $shapeC*$cosTheta;
      $poly[3] = $poly[1] - $shapeB*$sinTheta - $temp;
      $poly[9] = $poly[3] + 2*$temp;
      $poly[4] = $poly[2]*$fracHeight + $vertX*(1.0-$fracHeight);
      $poly[5] = $poly[3]*$fracHeight + $vertY*(1.0-$fracHeight);
      $poly[6] = $poly[8]*$fracHeight + $vertX*(1.0-$fracHeight);
      $poly[7] = $poly[9]*$fracHeight + $vertY*(1.0-$fracHeight);
      $point[-1]->[0] = $poly[0] - $backup*$cosTheta;
      $point[-1]->[1] = $poly[1] - $backup*$sinTheta;

      # round up the coords and indicate that the values are in pixels
      @poly = map { int($_) . 'p' } @poly;

      # draw the arrow
      $self->gd_polygon('filled', -coords => \@poly, -color => $color);
    }
  }

  foreach my $i ( 0 .. $#point - 1 ) {
    my ($sinTheta, $cosTheta);
    my $dx = $point[$i+1]->[0] - $point[$i]->[0];
    my $dy = $point[$i+1]->[1] - $point[$i]->[1];
    my $length = sqrt($dx*$dx + $dy*$dy);
    if ($length == 0) {
      $sinTheta = $cosTheta = 0.0;
    } else {
      $sinTheta = $dy/$length;
      $cosTheta = $dx/$length;
    }
    my $hw = $width/2.0;
    my $hwx = $hw*$sinTheta;
    my $hwy = $hw*$cosTheta;

    #
    #      1      dx    |
    #    $i   *        Theta
    #   4        *     /
    #      *        *
    #         *        *
    #   dy       *        2
    #               *   $i+1
    #                  3

    my @poly = (
		$point[$i]->[0] + $hwx, # 1
		$point[$i]->[1] - $hwy,
		$point[$i+1]->[0] + $hwx, # 2
		$point[$i+1]->[1] - $hwy,
		$point[$i+1]->[0] - $hwx, # 3
		$point[$i+1]->[1] + $hwy,
		$point[$i]->[0] - $hwx, # 4
		$point[$i]->[1] + $hwy,
	       );
    $self->image->gd_image->line(@{$point[$i]}, @{$point[$i+1]}, $self->image->color($color));
    # round up the coords and indicate that the values are in pixels
    @poly = map { int($_) . 'p' } @poly;

    # draw the line segment
    $self->gd_polygon('filled', -coords => \@poly, -color => $color);
  }
}

sub AUTOLOAD {
    my ($self, $arg) = (@_);
    $AUTOLOAD =~ /^.*::(\w+)$/; # strip the class name off the method name
    if ( $arg ) {
	$self->{"-$1"} = $arg;
	return $self;
    }
    else {
	return $self->{"-$1"};
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
package Image;
use GD;
use Data::Dumper;
use vars qw/@ISA %RGB/;
@ISA = qw/OGD::Object/;

{
  # Encapsulated class data
  my %_rgb;
  my %_color;

  my $file = "/usr/X11R6/lib/X11/rgb.txt";
  if ( -f $file ) {
    open RGB, "<$file";
    while (<RGB>) {
      chomp;
      if ( /^\s*(\d+)\s+(\d+)\s+(\d+)\s+(.+)$/ ) {
	my ($r, $g, $b, $key) = ($1, $2, $3, $4);
	$key =~ s/ //g;
	$key =~ tr/A-Z/a-z/;
	$_rgb{$key} = [$r, $g, $b];
      }
    }
    close RGB;
  }

  # Class methods, to operate on encapsulated class data
  sub _rgb {
    my ($self, $key) = @_;
    $key =~ s/ //g;
    $key =~ tr/A-Z/a-z/;
    exists $_rgb{$key} ? @{$_rgb{$key}} : ()
  }

  sub _color {
    my ($self, $key, $value) = @_;
    $key or die "undefined color key";
    if ( defined $value ) { # $value is an integer GD color index; can be zero
      $_color{$key} = $value;
      return $value;
    }
    else {
      return $_color{$key};
    }
  }
}

sub new {
  my ($class, $w, $h) = @_;
  my $self = {};
  bless $self, $class;
  $self->gd_image( new GD::Image($w, $h) );
  $self->w($w);
  $self->h($h);
  $self->gd_image->transparent($self->color('white'));
  return $self;
}

sub Frame {
  my ($self) = @_;
  return Frame->new(-parent => $self, -coords => [0, 0, $self->w - 1, $self->h - 1]);
}

sub color {
  my ( $self, $color ) = @_;
  my ($r, $g, $b);
  if ( ref($color) eq 'ARRAY' ) {
    ($r, $g, $b) = map {int $_} @$color;
  }
  elsif ($self->_rgb($color)) {
    ($r, $g, $b) = $self->_rgb($color);
  }
  elsif ( $color =~ /\#?([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})/ ) {
    ($r, $g, $b) = (hex($1), hex($2), hex($3));
  }
  else {
    die "unknown color: '$color'";
  }

  if ( $self->_color("$r.$g.$b") ) {
    return $self->_color("$r.$g.$b");
  }
  else {
    return $self->_color("$r.$g.$b", $self->gd_image->colorAllocate($r, $g, $b));
  }
}

sub hsv_to_rgb {
    # This procedure converts an HSB value to RGB.  It takes hue,
    # saturation, and value components (floating-point, 0-1.0) as
    # arguments, and returns a list containing RGB components
    # (integers, 0-65535) as result.  The code here is a copy of the
    # code on page 616 of "Fundamentals of Interactive Computer
    # Graphics" by Foley and Van Dam.

    my($self, $hue, $sat, $value) = @_;

    my($v, $i, $f, $p, $q, $t);

    $v = int(65535 * $value);
    return ($v, $v, $v) if $sat == 0;
    $hue *= 6;
    $hue = 0 if $hue >= 6;
    $i = int($hue);
    $f = $hue - $i;
    $p = int(65535 * $value * (1 - $sat));
    $q = int(65535 * $value * (1 - ($sat * $f)));
    $t = int(65535 * $value * (1 - ($sat * (1 - $f))));
    return ($v, $t, $p) if $i == 0;
    return ($q, $v, $p) if $i == 1;
    return ($p, $v, $t) if $i == 2;
    return ($p, $q, $v) if $i == 3;
    return ($t, $p, $v) if $i == 4;
    return ($v, $p, $q) if $i == 5;

} # end hsvToRgb

sub fill {
  my ($self, $color) = @_;
  $self->gd_image->filledRectangle(0, 0, $self->w - 1, $self->h - 1, $self->color($color));
  $self;
}

sub gif {
  shift->gd_image->gif;
}

sub png {
  return shift->gif;
  #shift->gd_image->png;
}

sub jpeg {
  shift->gd_image->jpeg;
}

sub cgidata {
  my ($self, $type) = @_;
  $type ||= ($self->can('gif') ? 'gif' : undef) || ($self->can('png') ? 'png' : undef);
  $type =~ tr/A-Z/a-z/;
  "Content-type: image/$type\n\n" . $self->gd_image->$type;
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
package Frame;
use GD;
use Data::Dumper;
use vars '@ISA';
@ISA = qw/OGD::Object/;

sub new {
  my ($class, %args) = @_;
  my $self = {};
  bless $self, $class;
  $self->parent($args{-parent});
  if ( ref($self->parent) eq 'Image' ) {
    $self->image($self->parent);
  }
  if ( ref($self->parent) eq 'Frame' ) {
    $self->image($self->parent->image);
  }
  $self->coords($args{-coords});
  return $self;
}

sub Frame {
  my ($self) = @_;
  return Frame->new(-parent => $self, -coords => [$self->coords]);
}

sub coords {
  my ( $self, @arg ) = @_;
  if ( @arg ) {
    if ( ref $arg[0] ) {
      $self->{-coords} = [map {int $_} @{$arg[0]}];
    }
    else {
      $self->{-coords} = [map {int $_} @arg];
    }
  }
  else {
    return @{$self->{-coords}};
  }
}


sub place {
  my ($self, %arg) = @_;

  my @parent_coord;
  if ( ref $self->parent eq 'Image' ) {
    @parent_coord = (0, 0, $self->image->w - 1, $self->image->h - 1);
  }
  elsif ( ref $self->parent eq 'Frame' ) {
    @parent_coord = $self->parent->coords;
  }

  my @c = @parent_coord;
  if ( $arg{-coords} ) {
    @c =  $self->translate_coords(@{$arg{-coords}});
  }

  $self->coords( @c );
  return $self;
}

sub stroke {
  my ($self, $color) = @_;
  $self->rectangle(0, 0, 1, 1, $color);
  $self;
}

sub fill {
  my ($self, $color) = @_;
  $self->filledRectangle(0, 0, 1, 1, $color);
  $self;
}

sub gd_rectangle {
  my ($self, $x1, $y1, $x2, $y2, $color, $filled) = @_;
  my $rectangle_method = $filled ? 'filledRectangle' : 'rectangle';
  my @c = $self->coords;
  $x1 = $c[0] + $x1*($c[2] - $c[0]);
  $y1 = $c[1] + $y1*($c[3] - $c[1]);
  $x2 = $c[0] + $x2*($c[2] - $c[0]);
  $y2 = $c[1] + $y2*($c[3] - $c[1]);
  $self->image->gd_image->$rectangle_method($x1, $y1, $x2, $y2, $self->image->color($color));
}

sub rectangle {
  my ($self, $x1, $y1, $x2, $y2, $color) = @_;
  $self->gd_rectangle($x1, $y1, $x2, $y2, $color, undef);
}

sub filledRectangle {
  my ($self, $x1, $y1, $x2, $y2, $color) = @_;
  $self->gd_rectangle($x1, $y1, $x2, $y2, $color, 'filled');
}

sub gd_polygon {
  my ($self, $filled, %arg) = @_;
  my $polygon_method = $filled ? 'filledPolygon' : 'polygon';

  my $poly = new GD::Polygon;
  my @c = $self->translate_coords(@{$arg{-coords}});
  foreach my $i ( 0 .. @c/2 - 1 ) {
    $poly->addPt(@c[2*$i, 2*$i+1]);
  }
  $self->image->gd_image->$polygon_method($poly, $self->image->color($arg{-color}));
}

sub polygon {
  my ($self, %arg) = @_;
  $self->gd_polygon(undef, %arg);
}

sub filledPolygon {
  my ($self, %arg) = @_;
  $self->gd_polygon('filled', %arg);
}

sub gd_arc {
  my ($self, $cx, $cy, $w, $h, $start, $end, $color) = @_;
  my @c = $self->coords;
  $cx = $c[0] + $cx*($c[2] - $c[0]);
  $cy = $c[1] + $cy*($c[3] - $c[1]);
  $self->image->gd_image->arc($cx, $cy, $h, $w, $start, $end, $self->image->color($color));
}

sub circle {
  my ($self, %arg) = @_;
  my ($cx, $cy) = @{$arg{-coords}};
  my $diameter = $arg{-diameter};
  my $outline = $arg{-outline};
  $self->gd_arc($cx, $cy, $diameter, $diameter, 0, 360, $outline);
}


sub arrow {
  my ($self, %arg) = @_;
  my $at = $arg{-at} * 0.99;
  my $beg = $arg{-beg};
  my $end = $arg{-end};
  my $fill = $arg{-fill} || '000000';
  my $width = ($arg{-width} || 0)/2;

  # determine the pixel size
  my @c = $self->coords;
  my $width_in_pixels = $width ? ($c[3] - $c[1])/$width : 1;
  my $flare = 40 * $width / $width_in_pixels;
  my $sweep = 25 * $width / $width_in_pixels;

  my $coords;
  if ( $end > $beg ) {
    $coords = [
	       $beg, $at - $width,
	       $end - $sweep, $at - $width,
	       $end - $sweep, $at - $width - $flare,
	       $end, $at,
	       $end - $sweep, $at + $width + $flare,
	       $end - $sweep, $at + $width,
	       $beg, $at + $width,
	      ]
  }
  else {
    $coords = [
	       $beg, $at - $width,
	       $end + $sweep, $at - $width,
	       $end + $sweep, $at - $width - $flare,
	       $end, $at,
	       $end + $sweep, $at + $width + $flare,
	       $end + $sweep, $at + $width,
	       $beg, $at + $width,
	      ]
  }
  $self->filledPolygon( -coords => $coords, -color => $fill );
}

sub directed_segment {
  my ($self, %arg) = @_;
  my $at = $arg{-at};
  my $beg = $arg{-beg};
  my $end = $arg{-end};
  my $color = $arg{-color} || '000000';

  # determine the pixel size
  my @c = $self->coords;
  my $frame_height_in_pixels = ($c[3] - $c[1]);
  my $flare = 5 / ($c[3] - $c[1]);
  my $sweep = 9 / ($c[2] - $c[0]);

  my $coords;
  if ( $end > $beg ) {
    $coords = [
	       $beg, $at - $flare,
	       $beg, $at + $flare,
	       $beg, $at,
	       $end, $at,
	       $end - $sweep, $at - $flare,
	       $end, $at,
	       $end - $sweep, $at + $flare
	      ]
  }
  else {
    $coords = [
	       $beg, $at - $flare,
	       $beg, $at + $flare,
	       $beg, $at,
	       $end, $at,
	       $end + $sweep, $at - $flare,
	       $end, $at,
	       $end + $sweep, $at + $flare
	      ]
  }
  $self->line( -coords => $coords, -color => $color );
}

sub string {
  my ($self, %arg) = @_;
  my $text = $arg{-text};
  my $color = $arg{-color} || '000000';
  my ($x, $y) = @{$arg{-coords}};
  my @c = $self->coords;
  if ($x =~ /^(-?\d+)p$/i) {
    $x = $c[0] + $1;
  }
  else {
    $x = $c[0] + $x*($c[2] - $c[0]);
  }
  if ($y =~ /^(-?\d+)p$/i) {
    $y = $c[1] + $1;
  }
  else {
    $y = $c[1] + $y*($c[3] - $c[1]);
  }
  $self->image->gd_image->string(gdSmallFont, $x, $y, $text, $self->image->color($color));  
}

sub ruler {
  my ($self, %arg) = @_;
  my @c = @{$arg{-coords}};
  my @limit = @{$arg{-limits}};
  my $major = $arg{-major} || 10;
  my $minor = $arg{-major} || 5;
  my $color = $arg{-color} || '000000';
  my $gd_color = $self->image->color($color);

  $self->image->gd_image->line(@c, $gd_color);

  my $numeric_step = int ($limit[1] - $limit[0]) / $major;
  my $major_step = ($c[2] - $c[0]) / $major;
  my $minor_step = $major_step / $minor;
  my $x = $c[0];
  my $y = $c[1];
  my $n = $limit[0];
  foreach my $i (0 .. $major) {
    $self->image->gd_image->line($x, $y-10, $x, $y, $gd_color);
    $self->string(-coords => [int($x - 45).'p', "-10p"], -text => int($n), -color => 'navy');
    last if $i == $major;
    foreach my $j (0 .. $minor - 1) {
      my $xminor = $x + $minor_step * $j;
      $self->image->gd_image->line($xminor, $y-5, $xminor, $y, $gd_color);
    }
    $x += $major_step;
    $n += $numeric_step;
  }
}

1;
