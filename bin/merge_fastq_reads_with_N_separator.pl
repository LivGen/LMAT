#!/usr/bin/perl -w

use FileHandle;
# Example:  Designed to merge two paired fastq files
# merge_sra_reads_with_N_separator.pl  $i.lib1.fastq $i.lib2.fastq $i.mergeReadN.fasta

my ($fastq_file1,$fastq_file2,$out);
if( $#ARGV == 0 ) {
   ($fastq_file1,$fastq_file2,$out) = split(/\s+/,$ARGV[0]);
} else {
   $fastq_file1=$ARGV[0];
   $fastq_file2=$ARGV[1];
   $out=$ARGV[2];
}

my $N="N";

open OUT,">$out" || die "failed to open [$out]";

my $fh1 = FileHandle::new();
my $fh2 = FileHandle::new();
$fh1->open("$fastq_file1") or die "Can't open [$fastq_file1]: $!\n";
$fh2->open("$fastq_file2") or die "Can't open [$fastq_file2]: $!\n";

my $sequence = "";
my $id="";

my $id_new="";
#while (my $line = <INDATA>) {
my $entry_cnt=0; 
for(;;) {
   my @entry1=getFastqEntry($fh1);
   last if( $#entry1 == -1 ); 
   my @entry2=getFastqEntry($fh2);
   my $id1=(split(/ /,$entry1[0]))[0];
   my $id2=(split(/ /,$entry1[0]))[0];
   if( $id1 eq $id2 ) {
      print OUT "$entry1[0]\n";
      print OUT "$entry1[1]N$entry2[1]\n";
      print OUT "$entry1[2]\n";
      print OUT "$entry1[3]#$entry2[3]\n";
   } else {
      print "error could not match read pairs\n";
   }
   $entry_cnt+=1; 
}
$fh1->close() or warn  "Can't close $fastq_file1: $!\n";
$fh2->close() or warn  "Can't close $fastq_file2: $!\n";

close OUT;

print "Completed merging $entry_cnt fastq entries\n";

sub getFastqEntry {
   my ($fh)= @_;
   my @result=();
   my $line="";
   if(defined( $line = <$fh> )) {
      chomp($line);
      push(@result,$line);
      defined( $line = <$fh> ) or die "readline failed for entry2: $!";
      chomp($line);
      push(@result,$line);
      defined( $line = <$fh> ) or die "readline failed for entry3: $!";
      chomp($line);
      push(@result,$line);
      defined( $line = <$fh> ) or die "readline failed for entry4: $!";
      chomp($line);
      push(@result,$line);
   }
   return @result;
}
