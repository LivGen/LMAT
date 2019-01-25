#!/usr/bin/perl
#
use strict;

my @private_ids=("USDA-LLNL","USAMRIID-LLNL", "SFG", "Calvin Keeler", "CDC", "Eureka", "Imigene", "Keim NAU", "LANL", "LLNL", "NBACC", "NBFAC", "PIADC", "USAMRIID", "USDA", "UTMB", "UVIC");

my %save;
while(my $line =<STDIN>) {
   chomp($line);
   ## sequences without GenBank ID need
   ## not be proprietary if they're assemblies/concatenations
   ## of draft contigs
   my $id="";
   if( $line =~ /Glued fragments/ && $line =~ /\((.*?), whole/ ) {
      $id=$1;
   }
   if( $id && $save{$id} ) {
      print "found dup: [$line] [$save{$id}]\n";
   } elsif( $id ) {
      $save{$id}=$line;
   }
}
