#!/usr/bin/perl

use strict;

my %g2t;
my $file="taxonomy/gi_taxid_nucl.dmp";
open(FILE,$file) || die "fail [$file]\n";
while(my $line = <FILE>) {
   chomp($line);
   my ($gi,$ti)=split(/\t/,$line);
   $g2t{$gi}=$ti;
}
close(FILE);

while(my $line = <STDIN> ) {
   if( $line =~ /gi\|(\d+)/) {
      my $mg=$1;
      if( !$g2t{$mg} ) {
         print "error but continue: [$mg]\n";
         next;
      }
      print "$mg\t$g2t{$mg}\n";
   }
}
