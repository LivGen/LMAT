#!/usr/bin/perl -w

use strict;

my %save;
while(my $line =<STDIN>) {
   chomp($line);
   my ($gi,$ti)=split(/\s+/,$line);   
   $save{$gi}=$ti;
}

open(FILE,$ARGV[0]) || die " faile\n";
while(my $line =<FILE>) {
   chomp($line);
   if( $line =~ /gi\|(\d+)/ ) {
      my $gi=$1;
      my $ti=$save{$gi};
      if( !$ti) {
         #print "no ti error $gi [$line]\n";
         next;
      }
      print "$ti\t$gi\t-1\t$line\n";
   } else {
      print "Error no gi: $line\n";
   }
}
