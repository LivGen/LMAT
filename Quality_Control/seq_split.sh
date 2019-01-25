#! /bin/sh 

#This little portion of code takes the list of seqid's being examined (dir.vir) and does the split of fragments with 200 bp and 50 overlap.    

for i in $(cat dir.ls);do
  cd "$i";
 /usr/gapps/kpath/bin/splitseq.pl "$i".fa "$i".fa.frag 200 50;
  cd ../
done  
 
    
