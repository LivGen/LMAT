#!/bin/sh

ramdisk=/l/ssd
cd $ramdisk
oname=10182017
size=$2
cutoff=$3

#Directory where taxonomy is:
#rdir=/p/lustre1/allen99/metagenomics/ref_sets/10182017update/taxonomy
rdir=./

#These two are needed
human_kmers=/p/lustre1/allen99/metagenomics/ref_sets/human-kmer-dumps/grand.kmerdump
adapt_kmers=/p/lustre1/allen99/metagenomics/ref_sets/artificial.seqs.dump
#Output directory
odir=$ramdisk
#hcnt=2512419327
hcnt=3000000000
nr=$rdir/numeric_rank.txt
/usr/gapps/kpath/lmat/LMAT-1.2.4b/bin/make_db_table -i /p/lustre2/bioinf/10182017update_files/$oname.flst -l -o $odir/$oname.20.db -k 20 -f $rdir/32To16.map -s $size -j $human_kmers -c $hcnt -u $adapt_kmers -g $cutoff -m $nr

#Move to directory where all the other databases are stored, if you have your own then change the path
cp $oname.20.db /p/lustre2/bioinf/dbs

cd $rdir
