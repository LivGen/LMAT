#!/bin/sh
threads=$1
oname=$2
devdir=/p/lscratchd/allen99/lmat-dev/preprocessing
seq 0 $threads | parallel $devdir/../rmodel/countTaxidFrequency -i ${oname}_th.{.} -f 32 -o $oname.kcnt.{.}
find . -maxdepth 1 -name $oname.kcnt.\*kcnt | python $devdir/../rmodel/combine_counts.py > $oname.kcnts
