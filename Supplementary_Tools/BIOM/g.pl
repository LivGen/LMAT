#!/usr/bin/perl -w

use strict;

## get the name of each file (line by line) from list_of_species_files
while(my $file=<STDIN>) {
    ### parse the file name to retrieve the RunID
    if( $file =~ /fc\.(.*?)_M_001\.fasta/ ) {
             my $fileid=$1;  ## this pulls out the SRR813269 in the filename             

             ## Now open up the file and print out the contents and add the RunID to each line            
            open(P,$file);
             <P>;
             # read each line from the file            
             while(my $res=<P>) { chomp $res;               
               ## print out RunID with the rest of the line of the file                 
               print "$fileid\t$res\n";              
             }              
             close(P);    
      }
   }


