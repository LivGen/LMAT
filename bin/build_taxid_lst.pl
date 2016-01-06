#!/usr/bin/perl -w 

use strict;

my %rank=readRank($ARGV[0]);
my $file_name=$ARGV[1];
my $ofile=$ARGV[2];
my $sstr=$ARGV[3];
procAll($file_name,\%rank,$ofile,$sstr);


sub procAll {
   my ($file,$rankRef,$ofile,$sstr) = @_;
   my %rank = %$rankRef;
   my %save_taxid;
   my %once;
   my (%cnt_call,%child,%parent);
   my $first=0;
   open(OFILE,">$ofile") || die "failed $ofile\n";
   open(FILE,$file) || die "failed $file\n";
   while(my $line = <FILE>) {
      chomp($line);   
      my @vals=split(/\t/,$line);
      my $taxid=$vals[2];
      if( !$rank{$taxid} ) {
         print "not in lineage file [$taxid]\n";
         next;
      }
      if( $rankRef->{$taxid} =~ /$sstr/ ) {
         print OFILE " " if( $first );
         print OFILE "$taxid";
         $first=1;
      }
   }
   close(FILE);
   print OFILE "\n" if( $first );
   close(OFILE);
   
}
   
sub readRank {
   my($file)=@_;
   my %res;
   open(FILE,$file) || die "fail $file\n";
   while(my $line =<FILE>) {
      chomp($line);
      $line =~ s/no rank/no_rank/g;
      my @vals=split(/\t/,$line);
      my ($ntid,$ktid)=("","");
      if( $vals[0]=~ /,taxid=(\d+),ktaxid=(\d+),/) {
         ($ntid,$ktid)=($1,$2);
      }
      $res{$ntid}=$line;
   }
   close(FILE);
   return %res;
}
