#!/usr/bin/env perl

# PODNAME: program_name.pl
# ABSTRACT: Description

use warnings;
use strict;
use Getopt::Long;
use autodie;
use Pod::Usage;

# get options
my %options;
get_and_check_options();

my $x = 0;
my $y = 0;
my $width = 5;
my $height = 20;

print q{<?xml version="1.0" encoding="UTF-8"?>}, "\n";
print qq{<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="$width" height="${height}pt" viewBox="0 0 $width $height" version="1.1">}, "\n";
print q{<g id="surface1">}, "\n";
#print qq{<rect x="0" y="0" width="$width" height="$height" style="fill:rgb(100%,100%,100%);fill-opacity:1;stroke:none;"/>}, "\n";

my @colours = (
    "#FFAE00", "#D99400", "#B37A00",
    "#8C6000", "#00FFBB", "#00D99F",
    "#00B383", "#008C67", "#FFEE00",
    "#C6B800", "#8C8300", "#00A2FF",
    "#007DC6", "#00598C", "#FF66BD",
    "#D9007B", "#B30065", "#8C004F"
);

my $stroke_info;
if( $options{lines} ){
   $stroke_info = 'stroke:black; stroke-width:0.5';
} else {
   $stroke_info = 'stroke:none';
}

if( defined $options{colours} ){
   @colours = @{$options{colours}};
}
foreach my $colour ( @colours ){
    for( my $i = 0; $i < 5; $i++ ){
        my $left_top = qq{$x $y};
        my $right_top = join(q{ }, $x + $width,  $y );
        my $right_bottom = join(q{ }, $x + $width, $y + $height );
        my $left_bottom = join(q{ }, $x, $y + $height );
        
        print join(q{ }, qq{<path style="$stroke_info ;fill-rule:nonzero;fill: $colour;fill-opacity:1;" d=},
            qq{"M $left_top L $right_top L $right_bottom L $left_bottom Z "/>}, "\n" );
        
        $x = $x + $width;
    }
}

print "</g>\n";
print "</svg>\n";


################################################################################
# SUBROUTINES
#
=func subroutine_name

  Usage       : subroutine_name( arguments )
  Purpose     : 
  Returns     : 
  Parameters  : 
  Throws      : 
  Comments    : None

=cut

=func get_and_check_options

  Usage       : get_and_check_options()
  Purpose     : parse the options supplied to the script using GetOpt::Long
  Returns     : None
  Parameters  : None
  Throws      : 
  Comments    : The default option are
                help which print a SYNOPSIS
                man which prints the full man page
                debug
                verbose

=cut

sub get_and_check_options {
    
    GetOptions(
        \%options,
        'lines+',
        'colours:s@',
        'help',
        'man',
        'debug+',
        'verbose',
    ) or pod2usage(2);
    
    # Documentation
    if( $options{help} ) {
        pod2usage( -verbose => 0, -exitval => 1, );
    }
    elsif( $options{man} ) {
        pod2usage( -verbose => 2 );
    }
    
    print "Settings:\n", map { join(' - ', $_, defined $options{$_} ? $options{$_} : 'off'),"\n" } sort keys %options if $options{verbose};
}

__END__

=pod

=head1 NAME

program_name.pl

=head1 DESCRIPTION

Description

=cut

=head1 SYNOPSIS

    program_name.pl [options] input file | STDIN
        --help                  print this help message
        --man                   print the manual page
        --debug                 print debugging information
        --verbose               turn on verbose output


=head1 ARGUMENTS

=over

arguments

=back

=head1 OPTIONS

**Same for optional arguments.

=over

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=head1 DEPENDENCIES

None

=head1 AUTHOR

=over 4

=item *

Richard White <richard.white@sanger.ac.uk>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2015 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut