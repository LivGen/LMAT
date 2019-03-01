#!/usr/bin/perl
#Verify that accession_taxid is made from assembly_summary_genbank.txt columns accession,taxid
use strict;

my %a2t;
my %a3t;
my $file="accession_to_taxid";
my $acc;
open(FILE,$file) || die "fail [$file]\n";
readline(FILE); #skip header file
while(my $line = <FILE>){
   chomp($line);
   my @t=split(/\t/,$line);
   my $ac=$t[0];
   my $ti=$t[2];
   my $tii=$t[1];	
   $a2t{$ti}=$ac;
   $a3t{$tii}=$ac;

   #print("$ac\t$ti\n"); 
}
close(FILE);

my $accession="";
while(my $line = <STDIN> ){
    if( $line =~ /\[tax_node_id\s(\d+)\]/){
        $acc="";
        my $tax=$1;
        if( !$a2t{$tax} ){

            $acc=$a3t{$tax};
        }else{
            $acc=$a2t{$tax};
	}
    print "$acc\t$tax\n";    
   }
}


