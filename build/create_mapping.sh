#!/bin/sh

header_file=$1 
devdir=/usr/gapps/kpath/lmat/utils/dbbuild

cat $header_file|$devdir/g2t.pl > seq2.mapping

