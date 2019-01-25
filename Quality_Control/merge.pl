#!/usr/bin/perl -w

use strict;

my $file=$ARGV[0];
my $ofile=$ARGV[1];
my %save;
while(my $line = <STDIN>) {
   chomp($line);
   my ($nm,$strn,$spec) = split(/\t/,$line);
   $save{$strn}=$spec;
}

open(OFILE,">$ofile") || die "fail $ofile\n";
open(FILE,$file) || die "fail $file\n";
while(my $line = <FILE>) {
   chomp($line);
   if( $line =~ /(\d+)$/ ) {
      my $spec=$1;
      my $conv=$save{$spec}; 
      if( !$conv ){
         print "unexpected skip: $line\n";
         next;
      } else {
         print OFILE "$line\t$conv\n";
      }
   } else {
      print "erorr: $line\n";
      exit(0);
   }
}
close(FILE);
close(OFILE);

