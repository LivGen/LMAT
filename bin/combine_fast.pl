#!/usr/bin/env perl 

use strict;

my %dict;
my %ds;
while(my $fn = <STDIN>) {
   chomp($fn);
   my @vals=split(/ /,$fn);
   foreach my $val (@vals) {
      if(!$val ) {
         last;
      }
      open(FILE,$val) || die "fail $val\n";
      while(my $line = <FILE> ) {
        chomp($line);
        my @parts = split(/;/,$line); 
        my $id=$parts[0];
         #print "[$parts[0]][$parts[1]][$parts[2]]\n";
        if( $dict{$id} ) {
            $ds{$id} += $parts[2];
            $dict{$id} += $parts[1];
        } else {
            $ds{$id} = $parts[2];
            $dict{$id} = $parts[1];
        }
      }
      close(FILE);
   }
}

foreach my $k (keys %ds) {
   my $pval=$k;
   $pval =~ s/\s+/\t/;
   print "$ds{$k}\t$dict{$k}\t$pval\n";
}
