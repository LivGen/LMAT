#!/bin/sh 

bdir=/usr/gapps/kpath/lmat/LMAT-1.2.4b/bin

$bdir/gen_rand_mod.sh --num_bases=1000000000 --tax_histo_cnt=/p/lscratchh/bioinf/db_bd2/2015update.kcnts --odir=/usr/mic/post1/metagenomics/ref_sets/fasta/01012015update/sinceMar2014 --db_file=/l/ssd/2015update.20.db --threads=40 --read_range=50:1000:25 --tax_dir=/usr/mic/post1/metagenomics/ref_sets/fasta/01012015update/sinceMar2014/taxonomy2

