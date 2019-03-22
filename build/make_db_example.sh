#!/bin/sh

bdir=/usr/gapps/kpath/lmat/LMAT-1.2.6/bin/

human_kmers=/usr/mic/post1/metagenomics/ref_sets/human-kmer-dumps/grand.kmerdump

output_db=$1

## list of files with output from kmerprefix
kmer_lst=$2  
kmer_size=20
memory=500
prune_size=200

$bdir/make_db_table -l -i $kmer_lst -o $output_db  -k $kmer_size -s $memory -j $human_kmers -f /usr/mic/post1/metagenomics/runtime_inputs/latest/m9.32To16.map -g $prune_size -m /
usr/mic/post1/metagenomics/runtime_inputs/latest/numeric_ranks -c 2524617332 -u /usr/mic/post1/metagenomics/ref_sets/artificial.seqs.dump


