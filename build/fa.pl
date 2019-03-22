#!/bin/perl
use strict;

#This script takes the full fasta per category, creates a directory with individual folders per genome

my $file = $ARGV[0];
my $dir = "$file.dir";

open(my $fh, '< :encoding(UTF-8)',$file) or die "Could not open file 'filename' $!";

if(-e $dir){
  chdir $dir; 
}else{mkdir $dir; chdir $dir;} 

my $count=0;
my $f="vnf"; open(my $ch, '>>',$f);
while( my $row = <$fh>)
{    
       $count++;
       chomp $row;  
       # if ($row =~ m/sequence_id\s(\d+)/){ the 2017update headers dont have sequence_id
	# so that it matches catalog in 2017update I just created a counter
       if($count != 0){	
  	   if(!-e $count){ print{$ch} "$count\n";      
           mkdir "$count", 0770; chdir "$count/";} else{chdir "$count";}   
           my $seq=<$fh>; 
           my $filename = "$count.fa";
        
           open(my $sf,'>',$filename)|| "file not created\n" ;
           print {$sf} $seq."\n"; close $sf;        
       } else{
           print "parse error, fix $row\n";

          } 
        chdir "../"; 
  }
close P;   
