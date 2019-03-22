#! /bin/sh -xvf 
#ifile='/usr/mic/post1/metagenomics/ref_sets/fasta/01012015update/sinceMar2014/2015_full.fasta'
ifile=$1
threads=40
oname=2015update
prefix=3
i=$2
#cd /p/lscratchh/bioinf/db_bd2/tst
cd ./
/usr/gapps/kpath/lmat/LMAT-1.2.4b/bin/kmerPrefixCounter -i $ifile.int -k 20 -o $oname -l $prefix  -f $i >& $oname.log$i



#seq 21 30| parallel /usr/gapps/kpath/lmat/LMAT-1.2.4b/bin/kmerPrefixCounter -i $ifile.int -k 20 -o $oname -l $prefix  -f {.} >& $oname.log
#sleep 3h

#seq 31 40| parallel /usr/gapps/kpath/lmat/LMAT-1.2.4b/bin/kmerPrefixCounter -i $ifile.int -k 20 -o $oname -l $prefix  -f {.} >& $oname.log
#sleep 3h



