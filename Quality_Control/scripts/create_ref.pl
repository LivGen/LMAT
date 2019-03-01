#! /bin/perl 

my $file = $ARGV[0];

open(FILE, $file) || die "Could not opem file\n";
print("TAXID\tsuperkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecie\n"); 
while(my $line=<FILE>)
{ my ($taxid,$kingdom,$superkingdom,$phylum, $class,$order,$family,$genus,$specie)=" ";
  chomp;
  my @t=split /\t/,$line;
  if($t[0]=~/,taxid=(\d+)/)
    {$taxid=$1;}
  if($line=~/\tsuperkingdom,(.*?)\t/){$superkingdom=$1;}elsif($line=~/\tsuperkingdom,(.*?)$/){$superkingdom=$1;}
  if($line=~/\tkingdom,(.*?)\t/){$kingdom=$1;}elsif($line=~/\tkingdom,(.*?)$/){$kingdom=$1;}
  if($line=~/\tphylum,(.*?)\t/){ $phylum=$1;}elsif($line=~/\tphylum,(.*?)$/){$phylum=$1;}
  if($line=~/\tclass,(.*?)\t/){ $class=$1;}elsif($line=~/\tclass,(.*?)$/){$class=$1;}
  if($line=~/\torder,(.*?)\t/){ $order=$1;}elsif($line=~/\torder,(.*?)$/){$order=$1;}
  if($line=~/\tfamily,(.*?)\t/){ $family=$1;}elsif($line=~/\tfamily,(.*?)$/){$family=$1;}
 if($line=~/\tgenus,(.*?)\t/){ $genus=$1;}elsif($line=~/\tgenus,(.*?)$/){$genus=$1;}
 if($line=~/\tspecies,(.*?)\t/){ $specie=$1;}elsif($line=~/\tspecies,(.*?)$/){$specie=$1;}
    
 # if($line=~/\tkingdom,(.*?)$/){$kingdom=$1;}
 # if($line=~/\tphylum,(.*?)$/){ $phylum=$1;}
 # if($line=~/\tclass,(.*?)$/){ $class=$1;}
 # if($line=~/\torder,(.*?)$/){ $order=$1;}
 # if($line=~/\tfamily,(.*?)$/){ $family=$1;}
# if($line=~/\tgenus,(.*?)$/){ $genus=$1;}
# if($line=~/\tspecies,(.*?)$/){ $specie=$1;}  
print("$taxid\t$superkingdom\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$specie\n");
}
close(FILE)

    
 

