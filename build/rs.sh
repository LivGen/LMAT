#!/bin/sh 

## example to run LMAT taxonomoy labeling 

## input can be fastq or fasta
file=$1
fname=`basename $file`

## assume output is in current directory
#$file.odir
dir=.  

bdir=/usr/gapps/kpath/lmat/LMAT-1.2.4b/bin

## command to put job in queue
if [ ! -e "$dir/$file.cleaned.db.lo.rl_output.0.30.fastsummary" ] ; then
   echo $fname
   sbatch --time=24:00:00 --di-mmap=npages=22000000,ver=stable -o  $dir/$fname.log.sbatch $bdir/catalyst_run_rl_ssd.sh --db_file=/dimmap/cleaned.db --query_file=$file --threads=96 --odir=$dir
fi
