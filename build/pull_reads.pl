#!/usr/bin/perl -w
use FileHandle;
use File::Basename;

my @myargs=split(/\s+/,$ARGV[0]);
my $ifile=$myargs[0];
chomp($ifile);
my $tbase=basename($ifile);
my $doShortRead=0;
my ($doLowScore,$min_score)=(0,0);
my $id_file=$myargs[1];
my $thresh=$myargs[2];
my $min_kmer=$myargs[3];
my $odir=$myargs[4];
my $select_type="exclusive"; ## make this the default
my $idname=basename($id_file);
my $ofilebase="$odir/$tbase.$idname.pulled";

sub writeFasta {
  my ($str,$header,$fh) = @_;
  print $fh ">$header\n";
  for(my $iter = 0; $iter < length($$str); $iter += 80) {
    my $mstr = substr($$str,$iter,80);
    print $fh "$mstr\n";
  }
}

if( !($select_type eq "exclusive") && !($select_type eq "inclusive") ){ 
   print "select type must be specified as exclusive or inclusive\n";
   exit(0);
}
my %ofhall;
my %valid;
open(FILE,$id_file) || die "fail to open id_file $id_file\n";
while(my $line=<FILE>) {
   chomp($line);
   if ( $select_type eq "inclusive" ) {
      $check{$line}=1;
   } else {
      my @vals=split(/\s+/,$line);
      if( $line =~ /^LowScore\s+/) {
         $doLowScore=1;
         $min_score=$vals[1];
         my $main_id="LowScore";
         my $ofile="$ofilebase.$main_id"; 
         $ofhall{$main_id}= FileHandle::new();
         $ofhall{$main_id}->open(">$ofile") || die "failed to create $ofile\n";
      } elsif( $line =~ /^ReadTooShort/) {
         $doShortRead=1;
         my $main_id="ReadTooShort";
         my $ofile="$ofilebase.$main_id"; 
         $ofhall{$main_id}= FileHandle::new();
         $ofhall{$main_id}->open(">$ofile") || die "failed to create $ofile\n";
      } else  {
         my $main_id=$vals[0];
         $valid{$main_id}=$main_id;
         for(my $it=1; $it <= $#vals; $it++) {
            $valid{$vals[$it]} = $main_id;
         }
         my $ofile="$ofilebase.$main_id"; 
         $ofhall{$main_id}= FileHandle::new();
         $ofhall{$main_id}->open(">$ofile") || die "failed to create $ofile\n";
      } 
   }
}
close(FILE);
my $fh;
if( $select_type eq "inclusive" ) {
   $fh = FileHandle::new();
   $fh->open(">$ofile") || die "fail to open ofile: $ofile\n";
}
my $cnt=0;
open(IFILE,$ifile) || die "failed to open $ifile\n";
while(my $line =<IFILE>) {
   chomp($line);
   my @vals=split(/\t/,$line);
   my ($tid,$score,$type)=split(/ /,$vals[4]);
   my ($ig1,$ig2,$valid_kmers) = split(/ /,$vals[2]);
   my ($ofh,$read)=("","");
   my $hdr="$vals[0];tid=$tid;score=$score;mtype=$type;valid_kmers=$valid_kmers;uid=$cnt;src=$tbase";
   if( $valid{$tid} && $score >= $thresh && $valid_kmers >= $min_kmer && !($vals[1] eq "X") )  {
     $read=$vals[1];
     $ofh = $ofhall{$valid{$tid}};
   } elsif( $doLowScore && $score < $min_score && $valid_kmers >= $min_kmer && !($vals[1] eq "X") ) {
     $read=$vals[1];
     $ofh = $ofhall{"LowScore"};
   } elsif( $type eq "NoDbHits" && $valid_kmers >= $min_kmer && !($vals[1] eq "X") ) {
     $read=$vals[1];
     $ofh = $ofhall{"NoDbHits"};
   } elsif( $type eq "ReadTooShort" && $doShortRead && !($vals[1] eq "X") ) {
     $read=$vals[1];
     $ofh = $ofhall{"ShortReads"};
   }
   if( $ofh ) {
      $cnt+=1;
      writeFasta(\$read,$hdr,$ofh);
   }
}
close(IFILE);

if( $select_type eq "exclusive" ) {
   foreach my $key (keys %check) {
      $check{$key}->close();
   }
} else {
   $fh->close();
}
