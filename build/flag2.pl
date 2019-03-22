#!/bin/perl -w
#input are the name of the dir to be evaluated and spec.cat* file
my $directory=$ARGV[0];
my $name;

##Retrieve genus information
if ($directory=~/all_(.*?).dir/){ $name=$1;}
my $gfile="genus.catalog.$name";
my $catalog="spec.catalog.$name"; 
my %gen1;


open(G, $gfile)|| die "Cant open $gfile";

while($_=<G>){ 
	chomp; my @s=split(/\t/,$_);
	my $seqid=$s[1]; my $genid=$s[3];
	if($genid=/,(\d+),/){
		$gen1{$seqid}=$1;

	}else{$gen1{$seqid}="-1";}

###  else everything regarding this genus is -1  
} close(G);  

open(my $sh,'<', $catalog)|| die "Couldnt open catalog";


#Open OutputFile
my $filename="summary_$name";
open(FILE1,'>>',$filename)||die "not created";
print FILE1 "Human_Frag\tTotIDFrag_Prop\tG_IDCounts%\tGFrag_Prop\t NoDb+Ls/frag\tIDPlasmid?\tSp_IDCounts-YES%\tSpFrag_Prop\tName\tSeqID\tTaxID\tSp_TaxID\tGen_TaxID\n"; #Print header for output file.

#Start computing percentages per $seqid
chdir $directory; 
while(my $catline=<$sh>){  
   
    chomp $catline;
    my @info= split /\t/, $catline;
    my $taxid=$info[3];
    if($info[3]=~/,(\d+),/){$taxid=$1;}
  #  print "$taxid\n";   
#Takes specie from spec.cat*
    my $seqid=$info[1];
   # print "$seqid\n";

###FILES IN USE
    my $basefile="$seqid/data.$seqid.fa.frag/$seqid.fa.frag.cleaned.db.lo.rl_output.0.30";
    my $fast="$basefile.fastsummary";
    my $fspecie="$fast.species";
    my $genfile="$fast.genus";
    my $pfile= "$fast.plasmid";
    my $ffrag="$seqid/$seqid.fa.frag";
    my $nms="$basefile.nomatchsum";   

###open file spec*-pct to attach information 
my $matched =0;my $count=0;my $totcount=0;my $line_cnt=0;my $hrc=0; my $frag=0; my $sfrag=0;

####GET NUMBER OF FRAGMENTS IN DOCUMENT###
open(F,$ffrag)||die "Could not open $ffrag";
while($_=<F>){
  if($_=~/>seq.(\d+)/){$frag++}else{next;}}
close(F);   
###CALCULATE SPECIES PERCENTAGE AND ETC 

if(-e $fspecie){
     #open *.SPECIE file and skip the header line
     open(S,'<',$fspecie)|| die "specie not found";
     <S>;
   
     while(my $row=<S>) 
     {  #If is not the header evaluate and store counts
        if ($row !~/^Average/){
         $line_cnt++;#We deal with the empty file
         chomp $row;
         my @comp = split(/\t/, $row);
	   
           	if ($taxid == $comp[3]){ #evaluates LMAT taxid  with specie call id
              	$count=$comp[2]; #stores correctly id counts
              	$totcount+=$comp[2]; }
                elsif($comp[3] eq 9606){ $hrc+=$comp[2];}
                elsif($taxid!=$comp[3] and $comp[3]!=32630){$totcount+=$comp[2];}
      #Calculate percentage
}}close(S);

         if($totcount!=0){ $matched=sprintf("%.2f",($count/$totcount)*100); $sfrag=sprintf("%.3f",($count/$frag));}
         else{$matched='NE';}

}else{$matched='NE';}
    
   
###OBSERVE IF IDENTIFIED AS PLASMID	
 my $st="FALSE";
if(-e $pfile ){
      
        open(P, $pfile)||print "Plasmid file not found";
        <P>; #skip header
        while(<P>){ my @comp=split(/\t/);
                   if($comp[3]==$taxid){$st="TRUE";}   
           } close P;
}else{ $st = "NE"};#files not exist   
          
###OPEN NODBHITS#####
my $ndb=0; my $ls=0, my $sc=0;
if(-e $nms){ open(N,$nms)||print "NoDBHits File not open: $nms";    
   while(my $info=<N>){ if($info=~/NoDbhits\t(\d+)/){$ndb=$1;}elsif($info=~/LowScore\t(\d+)/){$ls=$1;}}
   $sc=sprintf("%.2f",(($ndb+$ls)/$frag)*100);
} close(N);  

## Open genus information ###
my $gmatch=0; my $rc=0; my $gfrag=0;my $gtotrc=0;

if( -e $genfile){
open(GEN, $genfile)||die "Could not open $genfile";
<GEN>; #Skip header
while($_=<GEN>){ chomp; 
           my @t = split(/\t/, $_);
           if ($gen1{$seqid}==$t[3]){ #evaluates LMAT taxid  with gen call id
                $rc=$t[2]; #stores correctly id counts
                $gtotrc+=$t[2];}
           elsif($t[3]!=32630 and $t[3]!=$gen1{$seqid}){$gtotrc+=$t[2];} 
	   elsif($t[3]==9605){$hrc+=$t[2];}
            
      }close(GEN);

if($gtotrc !=0){
$gmatch=sprintf("%.2f",(($rc/$gtotrc)*100));
$gfrag=sprintf("%.3f",($rc/$frag))}; #Portion of fragments id at genus level 
}else{ $gmatch='NE'; $gfrag='NE';} 
 
#All Reads ID as something

my $ar=0;
if($fast){
open(Fo, $fast);
<Fo>;
while(<Fo>){chomp; my @l=split(/\t/); $ar+=$l[1];} close(Fo);
}else{ print "$seqid\n";}
 
#PRINT ALL INFORMATION

$ar=sprintf("%.2f",$ar/$frag); 
my $hfrag=sprintf("%.2f",($hrc/$frag)*100);
            
print FILE1 "$hfrag\t$ar\t$gmatch\t$gfrag\t$sc\t$st\t$matched\t$sfrag\t$catline\t$gen1{$seqid}\n"; 
}close(FILE1); close $sh;


