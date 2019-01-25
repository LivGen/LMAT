#!/bin/sh -vxf

nt_db=/p/lustre2/torres49/blastdb/db/nt
/usr/gapps/kpath/ncbi-blast-2.2.27+/bin/blastn -db $nt_db -query $1 -outfmt 7 -word_size 16 > ${1}.txt


