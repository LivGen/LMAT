#!/bin/perl

#This script checks if the fastsummary for the genome was done and if not it resubmits the job.

my $directory=$ARGV[0];
my $list="run.ls";
open(my $sh,'<', $list)|| die "Couldn't open list";
chdir $directory;
while(my $seqid=<$sh>){  
   
    chomp $seqid;
    my $fast="$seqid/data.$seqid.fa.frag/$seqid.fa.frag.cleaned.db.lo.rl_output.0.30.fastsummary";
    if(! -f $fast){
    	print "$seqid\n";
	chdir $seqid;
	system("/g/g19/bioinf/rs_dla.sh $seqid.fa.frag");
	chdir "../";
	}
 } close $sh;
