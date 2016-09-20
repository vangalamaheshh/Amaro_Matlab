#!/usr/bin/perl
use strict;
use warnings;

# notice to insert
my $notice = <<EON;
% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.
EON
#'

# main source files, destination directory
my @files = </xchip/gistic/GISTIC2.0/release_candidate_7/source/*.m>;
my $release_dir = '/xchip/gistic/GISTIC2.0/final_release/source/';

&notifiles;

# SegArray class source files, destination directory
@files = </xchip/gistic/GISTIC2.0/release_candidate_7/source/\@SegArray/*.m>;
$release_dir = '/xchip/gistic/GISTIC2.0/final_release/source/@SegArray/';

&notifiles;

exit 0;
###########################
sub notifiles {
    
    # loop across all m-files
    foreach( @files ) {
	if( m/.*\.m/ ) {
	    # process m-file 
	    my $file = $_;
	    print STDERR "$file\n";
	    # read file contents into memory
	    open F, "<$file" or die "couldn't open '$file' for reading.\n";
	    my @lines = <F>;
	    close F;
	    # remove path from file
	    $file =~ s/\/.*\///;
	    # initialize state machine
	    my $state = 0; #  0 = start;
	                   #  1 = first nonempty, noncomment ('live') line;
	                   #  2 = in comment
	                   #  3 = done (notice emitted)
	                   # write file line-by-line through state machine
	    open F, ">$release_dir$file" or die "couldn't open '$file' for writing.\n";
	    foreach( @lines ) {
		if( $state != 3 ) {
#	            print STDERR "$state ";
		    if( $state == 0 ) {
			if( m/^\s*%/ ) { # comment
			    $state = 2;
			}
			elsif( !m/^[ \t]*$/ and !m/\.\.\.\s*$/ ) {
                            # live uncontinued line
			    $state = 1;
			}
		    }
		    elsif( $state == 1 ) {
			if( m/^\s*%/ ) { # comment
			    $state = 2;
			}
			elsif( !m/^\s*$/ ) { # live line
			    # emit notice and exit state machine
			    print F "\n$notice\n";
			    $state = 3;
			}
		    }
		    elsif( $state == 2 ) {
			if( !m/^\s*%/ ) { # non-comment
			    # emit notice and exit state machine
			    if( m/^\s*$/ ) {
				print F "\n$notice\n";
			    }
			    else {
				# (supply extra blank line)
				print F "\n$notice\n";
			    }
			    $state = 3;
			}
		    }
		}
		# print source line 
		print F "$_";
	    }
	}
    } 
}
