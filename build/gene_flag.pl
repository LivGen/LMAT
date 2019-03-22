#!/bin/perl 
use strict;

# This script counts the amount of genes being reported for this genome and counts how many were actually labelled correctly to its specie.



#####GENESUMMARY SUMMARY##########

####METHOD TO GET UNIQUE GENES
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


#INPUT INFORMATIONN 
my $info=$ARGV[0];

my $dir="all_$info.dir";
my $fileS="summary_$info";
my $db="allgenes.7-14.20.db";

#Creates array to store the genes names
print "Seqid\tTotal_Count\tCountS&P\tprop\n";
open(F, $fileS)||die "Could not open summary file";
chdir "$dir";
<F>;
while(<F>){
my @t=split(/\t/,$_); 
my $seqid=$t[9];
my $spid=$t[11];#its envolved in commas  
if($spid=~/,(.*?),/){$spid=$1;}


my $genid=$t[12];
#open and read genesummary

my $fileG="$seqid/data.$seqid.fa.frag/tlst.$db.$seqid.fa.frag.$db.gl_output.0.1.20.genesummary";
if(-s $fileG and -e _){ 
my $tcnt=0; my $spcnt=0; my $pct=0;
my @genes;
my @spgenes;
my @filtered_genes; 

open(F2,$fileG)||die "Could not open $fileG\n";<F2>;
 while(<F2>){
  my @h=split(/\t/,$_);
#Puts sp level genes and gene level genes inside array
if($h[3]==$genid or $h[3]==$spid){
push(@spgenes,"$h[5]");
} #Fix to put name of gene inside the array

push(@genes, "$h[5]"); #insert all gene names 

#Count # of total non repetitive genes and total gene sp do fractional. 
}

#Filter unique genes
@spgenes=uniq(@spgenes);
@filtered_genes=uniq(@genes);
$tcnt=@filtered_genes;
$spcnt=@spgenes;

if($tcnt== 0){ $pct="NA";}else{ $pct=sprintf("%.3f",$spcnt/$tcnt);}
print("$seqid\t$tcnt\t$spcnt\t$pct\n");
}
}


