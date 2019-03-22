#!/usr/bin/perl -w

## Pull FASTQ reads for the headers specified in $headers

use strict;
my %sh;
my $headers=$ARGV[0];
open(FILE,$headers) || die "fail $headers\n";
while(my $line = <FILE>) {
   chomp($line);
   $sh{$line}=1;
}
close(FILE);

my $doPrn=0;
my $fq_file=$ARGV[1];
my $ofile=$ARGV[2];
open(OFILE,">$ofile") || die "fail $ofile\n";
open(FILE,$fq_file) || die "fail $fq_file\n";
while(my $line = <FILE>) {
   chomp($line);
   if( $line =~ /^@(.*)?/ ) {
      my $hid=$1;
      if( $sh{$hid} ) {
         $doPrn=1;
         print OFILE "$line\n";
      } else {
         $doPrn=0;
      }
   } elsif( $doPrn ) {
      print OFILE "$line\n";
   }
}
close(FILE);
close(OFILE);
