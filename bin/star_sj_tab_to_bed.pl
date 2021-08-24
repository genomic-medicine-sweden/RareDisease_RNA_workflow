#!/usr/bin/env perl

use warnings;
use strict;

my @motifs = qw{ non-canonical GT/AG CT/AC GC/AG CT/GC AT/AC GT/AT };
my @strands = qw{ . + - };


while(<>) {

    chomp;

    my @line_elements = split q{\t};

    next if $line_elements[6] == 0;
    
    ## Add contig, start position and end position 
    my @bed_elements = ( $line_elements[0], $line_elements[1] -1, $line_elements[2] );

    ## build name column
    my @names = ( 
	qq{motif=$motifs[$line_elements[4]]}, 
    	qq{uniquely_mapped=$line_elements[6]}, 
	qq{multi_mapped=$line_elements[7]},
	qq{maximum_spliced_alignment_overhang=$line_elements[8]} 
    );

    ## join and add name columns
    push @bed_elements, join q{;}, @names;

    ## Add score and strand
    push @bed_elements, $line_elements[6], $strands[$line_elements[3]];

    print STDOUT join(qq{\t}, @bed_elements) . qq{\n};
}
