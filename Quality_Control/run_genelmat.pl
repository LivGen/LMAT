#!/usr/bin/perl -w

use strict;
use FileHandle;
use File::Basename;
my $fh;
my $cnt=0;
#my $db="lmat-4-14.20mer.db";
## default local LMAT taxonomy DB name
my $db="cleaned.db";

## gene dbs:
## AMR predicted db
## amr.08252014.16.db
## AMR curated db
## 09092014.amr_cur.db
## default "ALL gene library"
#  allgenes.7-14.20.db
my $gdb="allgenes.7-14.20.db";

while(my $file=<STDIN>) {
   chomp($file) if($file);
   my $dir="\.\/$file\/data.$file.fa.frag";
    
  my $pip = "find `$dir` -name $file.fa.frag.$db.lo.rl_output\*out | ";
  open(PIP,$pip) || die "fail [$pip]\n";
  chdir $dir;
 my $fname="tlst.$gdb.$file.fa.frag";
  #open(F, $fname);  
#system("ls -1 \*out > tlst.$gdb.$file.fa.frag");
   open(OFILE,">$fname") ||die "fail [$fname]\n";
   while(my $lin=<PIP>) {
     print OFILE $lin;
   }
   close(OFILE);
   close(PIP);
  # close(F);
   if( -e "./$dir/$fname") {
      #my $cmd = "sbatch --time=24:00:00 --di-mmap=npages=20500000,ver=stable -o $fname.log /g/g21/allen99/lmatnotes/nasa/custom_catalyst_run_gl.sh --odir=$dir --db_file=/dimmap/$gdb --ilst=$fname --overwrite";
      my $cmd = "sbatch --time=24:00:00 --di-mmap=npages=20500000,ver=stable -o ./$fname.log catalyst_run_gl.sh --odir=$dir --db_file=/dimmap/$gdb --ilst=$fname --overwrite";
      print "$cmd\n";
     #system("$cmd");
     } else {
      print "Unexpected error $file\n";}
chdir "../../"; 
     
}



