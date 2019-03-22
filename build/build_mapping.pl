#!/usr/bin/perl
#
use strict;

my @private_ids=("USDA-LLNL","USAMRIID-LLNL", "SFG", "Calvin Keeler", "CDC", "Eureka", "Imigene", "Keim NAU", "LANL", "LLNL", "NBACC", "NBFAC", "PIADC", "USAMRIID", "USDA", "UTMB", "UVIC");

while(my $line =<STDIN>) {
   ## sequences without GenBank ID need
   ## not be proprietary if they're assemblies/concatenations
   ## of draft contigs
   my $private=0;
   my ($tid,$gid,$kid) = (-1,-1,-1);
   if( $line =~ /tax_node_id (\d+)/ ) {
      $tid=$1;
   }
   next if($tid >= 10000000);
   foreach my $id (@private_ids) {
      if( $line =~ / $id / ) {
         $private=1;
         last;
      }
   }
   if( $line =~ /gi\|(\d+)/ ) {
      $gid=$1;
   }
   if( $line =~ /sequence_id (\d+)/ ) {
      $kid=$1;
   }
   ## assume synthetic construct
   if( $tid == -1 ) {
      $tid=32630;
   }
   if( !$private ) {
      print "$tid\t$gid\t$kid\t$line";
   }
}
