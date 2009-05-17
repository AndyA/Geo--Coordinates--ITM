package Geo::Coordinates::ITM;

use warnings;
use strict;

use base qw( Exporter );

our @EXPORT_OK = qw( ll_to_grid grid_to_ll );

use Carp;
use Math::Trig;

=head1 NAME

Geo::Coordinates::ITM - Convert coordinates between lat/lon and Irish Transverse Mercator

=head1 VERSION

This document describes Geo::Coordinates::ITM version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

  use Geo::Coordinates::ITM;
  
=head1 DESCRIPTION

=head1 INTERFACE 

=head2 C<< ll_to_grid >>

=cut

sub ll_to_grid {
  my ( $lat, $long ) = @_;

  return (
    _ll2e(
      $lat,   $long,    6378137,  6356752.314,
      600000, 0.999820, 53.50000, -8.00000
    ),
    _ll2n(
      $lat,   $long,  6378137,  6356752.314,
      600000, 750000, 0.999820, 53.50000,
      -8.00000
    )
  );
}

=head2 C<< grid_to_ll >>

=cut

sub grid_to_ll {
  my ( $e, $n ) = @_;

  return (
    _en2lat(
      $e,     $n,     6378137,  6356752.314,
      600000, 750000, 0.999820, 53.50000,
      -8.00000
    ),
    _en2lon(
      $e,     $n,     6378137,  6356752.314,
      600000, 750000, 0.999820, 53.50000,
      -8.00000
    )
  );
}

sub _pow {
  my ( $a, $b ) = @_;
  return $a**$b;
}

sub _ll2e {
  my ( $PHI, $LAM, $a, $b, $e0, $f0, $PHI0, $LAM0 ) = @_;

  # Convert angle measures to radians
  my $Pi      = pi;
  my $RadPHI  = $PHI * ( $Pi / 180 );
  my $RadLAM  = $LAM * ( $Pi / 180 );
  my $RadPHI0 = $PHI0 * ( $Pi / 180 );
  my $RadLAM0 = $LAM0 * ( $Pi / 180 );
  my $af0     = $a * $f0;
  my $bf0     = $b * $f0;
  my $e2      = ( _pow( $af0, 2 ) - _pow( $bf0, 2 ) ) / _pow( $af0, 2 );
  my $n       = ( $af0 - $bf0 ) / ( $af0 + $bf0 );
  my $nu = $af0 / ( sqrt( 1 - ( $e2 * _pow( sin( $RadPHI ), 2 ) ) ) );
  my $rho = ( $nu * ( 1 - $e2 ) )
   / ( 1 - ( $e2 * _pow( sin( $RadPHI ), 2 ) ) );
  my $eta2 = ( $nu / $rho ) - 1;
  my $p    = $RadLAM - $RadLAM0;
  my $IV   = $nu * ( cos( $RadPHI ) );
  my $V
   = ( $nu / 6 ) 
   * ( _pow( cos( $RadPHI ), 3 ) )
   * ( ( $nu / $rho ) - ( _pow( tan( $RadPHI ), 2 ) ) );
  my $VI
   = ( $nu / 120 ) 
   * ( _pow( cos( $RadPHI ), 5 ) )
   * ( 5 
     - ( 18 * ( _pow( tan( $RadPHI ), 2 ) ) )
     + ( _pow( tan( $RadPHI ), 4 ) )
     + ( 14 * $eta2 )
     - ( 58 * ( _pow( tan( $RadPHI ), 2 ) ) * $eta2 ) );

  return $e0 + ( $p * $IV ) + ( _pow( $p, 3 ) * $V )
   + ( _pow( $p, 5 ) * $VI );
}

sub _ll2n {
  my ( $PHI, $LAM, $a, $b, $e0, $n0, $f0, $PHI0, $LAM0 ) = @_;

  my $Pi      = pi;
  my $RadPHI  = $PHI * ( $Pi / 180 );
  my $RadLAM  = $LAM * ( $Pi / 180 );
  my $RadPHI0 = $PHI0 * ( $Pi / 180 );
  my $RadLAM0 = $LAM0 * ( $Pi / 180 );
  my $af0     = $a * $f0;
  my $bf0     = $b * $f0;
  my $e2      = ( _pow( $af0, 2 ) - _pow( $bf0, 2 ) ) / _pow( $af0, 2 );
  my $n       = ( $af0 - $bf0 ) / ( $af0 + $bf0 );
  my $nu = $af0 / ( sqrt( 1 - ( $e2 * _pow( sin( $RadPHI ), 2 ) ) ) );
  my $rho = ( $nu * ( 1 - $e2 ) )
   / ( 1 - ( $e2 * _pow( sin( $RadPHI ), 2 ) ) );
  my $eta2 = ( $nu / $rho ) - 1;
  my $p    = $RadLAM - $RadLAM0;
  my $M    = _marc( $bf0, $n, $RadPHI0, $RadPHI );
  my $I    = $M + $n0;
  my $II   = ( $nu / 2 ) * ( sin( $RadPHI ) ) * ( cos( $RadPHI ) );
  my $III
   = (
    ( $nu / 24 ) * ( sin( $RadPHI ) ) * ( _pow( cos( $RadPHI ), 3 ) ) )
   * ( 5 - ( _pow( tan( $RadPHI ), 2 ) ) + ( 9 * $eta2 ) );
  my $IIIA
   = (
    ( $nu / 720 ) * ( sin( $RadPHI ) ) * ( _pow( cos( $RadPHI ), 5 ) ) )
   * ( 61 
     - ( 58 * ( _pow( tan( $RadPHI ), 2 ) ) )
     + ( _pow( tan( $RadPHI ), 4 ) ) );

  return $I + ( _pow( $p, 2 ) * $II ) + ( _pow( $p, 4 ) * $III )
   + ( _pow( $p, 6 ) * $IIIA );
}

sub _en2lat {
  my ( $East, $North, $a, $b, $e0, $n0, $f0, $PHI0, $LAM0 ) = @_;

  my $Pi      = pi;
  my $RadPHI0 = $PHI0 * ( $Pi / 180 );
  my $RadLAM0 = $LAM0 * ( $Pi / 180 );

  my $af0 = $a * $f0;
  my $bf0 = $b * $f0;
  my $e2  = ( _pow( $af0, 2 ) - _pow( $bf0, 2 ) ) / _pow( $af0, 2 );
  my $n   = ( $af0 - $bf0 ) / ( $af0 + $bf0 );
  my $Et  = $East - $e0;

  my $PHId = _init_lat( $North, $n0, $af0, $RadPHI0, $n, $bf0 );

  my $nu = $af0 / ( sqrt( 1 - ( $e2 * ( _pow( sin( $PHId ), 2 ) ) ) ) );
  my $rho
   = ( $nu * ( 1 - $e2 ) ) / ( 1 - ( $e2 * _pow( sin( $PHId ), 2 ) ) );
  my $eta2 = ( $nu / $rho ) - 1;

  my $VII = ( tan( $PHId ) ) / ( 2 * $rho * $nu );
  my $VIII
   = ( ( tan( $PHId ) ) / ( 24 * $rho * _pow( $nu, 3 ) ) )
   * ( 5 
     + ( 3 * ( _pow( tan( $PHId ), 2 ) ) ) 
     + $eta2
     - ( 9 * $eta2 * ( _pow( tan( $PHId ), 2 ) ) ) );
  my $IX
   = ( ( tan( $PHId ) ) / ( 720 * $rho * _pow( $nu, 5 ) ) )
   * ( 61 
     + ( 90 * ( ( tan( $PHId ) ) ^ 2 ) )
     + ( 45 * ( _pow( tan( $PHId ), 4 ) ) ) );

  my $_en2lat
   = ( 180 / $Pi )
   * ( $PHId 
     - ( _pow( $Et, 2 ) * $VII ) 
     + ( _pow( $Et, 4 ) * $VIII )
     - ( ( $Et ^ 6 ) * $IX ) );

  return ( $_en2lat );
}

sub _en2lon {
  my ( $East, $North, $a, $b, $e0, $n0, $f0, $PHI0, $LAM0 ) = @_;

  my $Pi      = pi;
  my $RadPHI0 = $PHI0 * ( $Pi / 180 );
  my $RadLAM0 = $LAM0 * ( $Pi / 180 );

  my $af0 = $a * $f0;
  my $bf0 = $b * $f0;
  my $e2  = ( _pow( $af0, 2 ) - _pow( $bf0, 2 ) ) / _pow( $af0, 2 );
  my $n   = ( $af0 - $bf0 ) / ( $af0 + $bf0 );
  my $Et  = $East - $e0;

  my $PHId = _init_lat( $North, $n0, $af0, $RadPHI0, $n, $bf0 );

  my $nu = $af0 / ( sqrt( 1 - ( $e2 * ( _pow( sin( $PHId ), 2 ) ) ) ) );
  my $rho
   = ( $nu * ( 1 - $e2 ) ) / ( 1 - ( $e2 * _pow( sin( $PHId ), 2 ) ) );
  my $eta2 = ( $nu / $rho ) - 1;

  my $X = ( _pow( cos( $PHId ), -1 ) ) / $nu;
  my $XI = ( ( _pow( cos( $PHId ), -1 ) ) / ( 6 * _pow( $nu, 3 ) ) )
   * ( ( $nu / $rho ) + ( 2 * ( _pow( tan( $PHId ), 2 ) ) ) );

  my $XII
   = ( ( _pow( cos( $PHId ), -1 ) ) / ( 120 * _pow( $nu, 5 ) ) )
   * ( 5 
     + ( 28 * ( _pow( tan( $PHId ), 2 ) ) )
     + ( 24 * ( _pow( tan( $PHId ), 4 ) ) ) );

  my $XIIA
   = ( ( _pow( cos( $PHId ), -1 ) ) / ( 5040 * _pow( $nu, 7 ) ) )
   * ( 61 
     + ( 662 *  ( _pow( tan( $PHId ), 2 ) ) )
     + ( 1320 * ( _pow( tan( $PHId ), 4 ) ) )
     + ( 720 *  ( _pow( tan( $PHId ), 6 ) ) ) );

  my $_en2lon
   = ( 180 / $Pi )
   * ( $RadLAM0 
     + ( $Et * $X ) 
     - ( _pow( $Et, 3 ) * $XI )
     + ( _pow( $Et, 5 ) * $XII ) 
     - ( _pow( $Et, 7 ) * $XIIA ) );

  return $_en2lon;
}

sub _init_lat {
  my ( $North, $n0, $afo, $PHI0, $n, $bfo ) = @_;

  my $PHI1 = ( ( $North - $n0 ) / $afo ) + $PHI0;

  my $M = _marc( $bfo, $n, $PHI0, $PHI1 );

  my $PHI2 = ( ( $North - $n0 - $M ) / $afo ) + $PHI1;

  while ( abs( $North - $n0 - $M ) > 0.00001 ) {
    $PHI2 = ( ( $North - $n0 - $M ) / $afo ) + $PHI1;
    $M = _marc( $bfo, $n, $PHI0, $PHI2 );
    $PHI1 = $PHI2;
  }

  return $PHI2;
}

sub _marc {
  my ( $bf0, $n, $PHI0, $PHI ) = @_;
  return $bf0 * (
    (
      (
           1 
         + $n 
         + ( ( 5 / 4 ) * _pow( $n, 2 ) )
         + ( ( 5 / 4 ) * _pow( $n, 3 ) )
      ) * ( $PHI - $PHI0 )
    ) - (
      (
          ( 3 * $n ) 
        + ( 3 * _pow( $n, 2 ) )
         + ( ( 21 / 8 ) * _pow( $n, 3 ) )
      ) * ( sin( $PHI - $PHI0 ) ) * ( cos( $PHI + $PHI0 ) )
     ) + (
      (
        ( ( 15 / 8 ) * _pow( $n, 2 ) ) + ( ( 15 / 8 ) * _pow( $n, 3 ) )
      ) 
      *  ( sin( 2 * ( $PHI - $PHI0 ) ) )
       * ( cos( 2 * ( $PHI + $PHI0 ) ) )
     ) - (
      ( ( 35 / 24 ) * _pow( $n, 3 ) ) 
      *  ( sin( 3 * ( $PHI - $PHI0 ) ) )
       * ( cos( 3 * ( $PHI + $PHI0 ) ) )
     )
  );
}

1;
__END__

=head1 BUGS

Please report any bugs or feature requests to
C<bug-geo-coordinates-itm@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.

=head1 AUTHOR

Andy Armstrong  C<< <andy@hexten.net> >>

Code gratefully stolen from L<http://bit.ly/ZqpUA>

  http://svn.geograph.org.uk/svn/branches/british-isles/libs/geograph \
    /conversionslatlong.class.php

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2009, Andy Armstrong C<< <andy@hexten.net> >>.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.
