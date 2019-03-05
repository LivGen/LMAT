#!/bin/sh -xvf

#This script assumes that the input file is the concatenation of all the *species files up for observation
#The first column should be the shortened filenames -> ls -1 *.species| ~/g.pl >> new.allspecies  

filename=$1

if [ -z $2 ];then
     threshold=0
else
     threshold=$2
fi 

/usr/bin/Rscript BIOM.R $filename $threshold --save
