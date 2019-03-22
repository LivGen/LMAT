#!/usr/bin/perl -w

use strict;

my $startid=9000000;
my %cnt;

while(my $line = <STDIN>) {
   my @vals=split(/\t/,$line);
   $cnt{$vals[0]} .= "$line";
}

my %humanIds;
foreach my $val (keys %cnt) {
   my @vals=split(/\n/,$cnt{$val});
   if( $#vals > 0 ) {
      foreach my $chk (@vals) {
         my $tid=(split(/\t/,$chk))[0];
         my $nid=$tid;
         $startid++;
         $nid=$startid;
         print "$tid\t$chk\n";
      }
   } else {
         my $tid=(split(/\t/,$vals[0]))[0];
         print "$tid\t$vals[0]\n";
   }
}
