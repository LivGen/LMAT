#!/bin/sh -xvf

#Remember to export the LMAT_DIR

start=275
end=1000
bdir=/usr/gapps/kpath/lmat/LMAT-1.2.4b/bin
wdir=/p/lustre1/allen99/metagenomics/runtime_inputs/2015update
num_bases=1000000000
tax_histo_cnt=$wdir/2015update.kcnts


for i in `seq $start 25 $end`; do
	sbatch --time=24:00:00 $bdir/gen_rand_mod.sh --num_bases=$num_bases --odir=/p/lustre2/bioinf/db_bd2 --db_file=$wdir/2015update.20.db --threads=40 --read_len=$i --tax_dir=$wdir/taxonomy ;
done 



