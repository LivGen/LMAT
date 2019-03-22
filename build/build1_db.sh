#!/bin/sh -xvf
ifile=$1
threads=$2
oname=$3
prefix=$4
## note, 16 bit is not yet support in the code
tid_bits=32
/usr/gapps/kpath/lmat/utils/dbbuild/build_header_table_v2.py $ifile seq.mapping.multi .
seq 0 $threads | parallel /usr/gapps/kpath/lmat/LMAT-1.2.4b/bin/kmerPrefixCounter -i $ifile.int -k 20 -o $oname -l $prefix -f {.} >& $oname.log
seq 0 $threads | parallel /usr/gapps/kpath/lmat/LMAT-1.2.4b/bin/tax_histo -f $tid_bits -o ${oname}_th.{.} -d $oname.{.} -t taxonomy/ncbi_taxonomy.segment.dat  >& $oname.log
find . -maxdepth 1 -name ${oname}_th\* | sort -t"." -k2n,2n > $oname.flst
