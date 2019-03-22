#!/usr/bin/perl -w 

use strict;

my @myargs=split(/\s+/,$ARGV[0]);
my $file_name=$myargs[0];
my $sig_thresh=$myargs[1];
my $min_kmers=$myargs[2];
my $min_tax_score=$myargs[3];

procAll($file_name,$sig_thresh,$min_kmers,$min_tax_score); 


sub procAll {
   my ($file,$thresh,$min_kmers,$min_tax_score) = @_;
   my %save_taxid;
   my %once;
   my $valid_read=1;
   my (%cnt_call,%cnt_call_sum);
   my (%tax_cnt_call,%tax_cnt_call_sum);
   open(FILE,$file) || die "failed $file\n";
   while(my $line = <FILE>) {
      chomp($line);
      my @vals=split(/\t/,$line);
      #there seems to be a tab that slips in here, must adjust
      my ($idx1,$idx2,$idx3)=(3,2,4);
      if( !$vals[3]) {
         ($idx1,$idx2,$idx3)=(3+1,2,4+1);
      }
      if( $min_kmers != -1 ) {
         if( $vals[$idx1] ) {
            my $valid_kmers = (split(/ /,$vals[$idx1]))[2];
            if( $valid_kmers < $min_kmers) {
               next; 
            }
         } else {
            print "Implies bad read [$vals[3]] [$line]\n";
            print "check [$vals[0]][$vals[1]][$vals[2]][$vals[3]][$vals[4]][$vals[5]]\n";
            exit(0);
         }
      }
      my ($ktaxid,$tax_score) = split(/ /,$vals[$idx2]);
      my ($gid,$gs,$ignore) = split(/ /,$vals[$idx3]);

      my $id="$ktaxid $gid";
      if( $gs >= $thresh ) {
         $cnt_call{$id}++;
         $cnt_call_sum{$id} += $gs;
         if( $tax_score >= $min_tax_score ) {
            $tax_cnt_call{$id}++;
            $tax_cnt_call_sum{$id} += $gs;
         }
      }
   }
   close(FILE);

   open(OFILE,">$file.$sig_thresh.$min_kmers.genesummary") || die "failed to write to $file.summary";
   foreach my $node (keys %cnt_call) {
      my $score = $cnt_call_sum{$node};
      my $cnt=$cnt_call{$node};
      my ($tid,$gid)=split(/ /,$node);
      print OFILE "$score\t$cnt\t$tid\t$gid\n";
   }
   close(OFILE);
   open(OFILE1,">$file.$sig_thresh.$min_kmers.genesummary.min_tax_score.$min_tax_score") || die "failed to write to $file.summary";
   foreach my $node (keys %tax_cnt_call) {
      my $score = $tax_cnt_call_sum{$node};
      my $cnt=$tax_cnt_call{$node};
      my ($tid,$gid)=split(/ /,$node);
      print OFILE1 "$score\t$cnt\t$tid\t$gid\n";
   }
   close(OFILE1);
}
