#!/usr/bin/perl -w 

use strict;

my @myargs=split(/\s+/,$ARGV[0]);
my $file_name=$myargs[0];
my %rank=readRank($myargs[1]);
my $sig_thresh=$myargs[2];
my $min_kmers=$myargs[3];

procAll($file_name,\%rank,$sig_thresh,$min_kmers); 


sub procAll {
   my ($file,$rankRef,$thresh,$min_kmers) = @_;
   my %rank = %$rankRef;
   my %save_taxid;
   my %once;
   my $valid_read=1;
   my (%cnt_call,%cnt_call_sum, %child,%parent);
   open(FILE,$file) || die "failed $file\n";
   while(my $line = <FILE>) {
      chomp($line);
      my @vals=split(/\t/,$line);
      if( $min_kmers != -1 ) {
         my $valid_kmers = (split(/ /,$vals[2]))[2];
         if( $valid_kmers < $min_kmers) {
            $cnt_call{"ShortRead"}++;
            next; 
         }
      }
      my $lidx=$#vals;
      my ($ktaxid,$loscore,$label_type) = split(/ /,$vals[4]);
      if( !$label_type ) {
         print "huh: $line\n";
         print "$vals[4]\n";
         print "$vals[3]\n";
         exit(0);
      }
      if( $label_type eq "ReadTooShort" ) {
         $cnt_call{"ShortRead"}++;
         next;
      }
      if($label_type eq "NoDbHits" ) {
         $cnt_call{"NoMatch"}++;
         next;
      }
      if($label_type eq "LCA_ERROR" ) {
         $cnt_call{"LCA_ERROR"}++;
         next;
      }
      if($loscore < $thresh ) {
         $cnt_call{"LowScore"}++;
         next;
      }
      if($ktaxid <= 0 ) {
         print "When taxid=[$ktaxid] is <= 0, a previous category should have been triggered\n";
         print "$line\n";
         exit(0);
      }
      my $lineage=$rank{$ktaxid};
      if( !$lineage) {
         if( !$once{$ktaxid} ) {	
            $once{$ktaxid} = 1; 
         }
         next;
      }
      my @lv=split(/\t/,$lineage);
      if( $#lv == 0 ) {
         my $curr="root"; 
         $cnt_call{$curr}++;
         $cnt_call_sum{$curr} += $loscore;
         $save_taxid{$curr} = $ktaxid;
      } else {
         my $curr=$lv[$#lv];
         $cnt_call{$curr}++;
         $cnt_call_sum{$curr} += $loscore;
         $save_taxid{$curr} = $ktaxid;
      }
   }
   close(FILE);

   open(OFILE,">$file.$sig_thresh.$min_kmers.fastsummary") || die "failed to write to $file.summary";
   open(OFILE1,">$file.$sig_thresh.$min_kmers.nomatchsum") || die "failed to write to $file.summary";
   foreach my $node (keys %cnt_call) {
      if( !$child{$node} ) {
         my $leaf_node_sum=$cnt_call_sum{$node};
         my $leaf_node_cnt=$cnt_call{$node};
         my $tot_read_cnt=$leaf_node_cnt;
         my $lstr="";
         if( !$save_taxid{$node} ) {
            $lstr="$node;$leaf_node_cnt";
            print OFILE1 "$lstr\n";
         } else {
            $lstr="$save_taxid{$node}\t$node;$leaf_node_cnt;$leaf_node_sum";
            print OFILE "$lstr\n";
         }


      }
   }
   close(OFILE);
   close(OFILE1);
}
   
sub readRank {
   my($file)=@_;
   my %res;
   open(FILE,$file) || die "fail $file\n";
   while(my $line =<FILE>) {
      chomp($line);
      $line =~ s/no rank/no_rank/g;
      my @vals=split(/\t/,$line);
      my $tid="";
      if( $vals[0]=~ /ktaxid=(\d+),/) {
         $tid=$1;
      } else {
         print "Error failed to parse taxonomy lineage: $line\n";
      }
      $res{$tid}=$line;
   }
   close(FILE);
   return %res;
}
