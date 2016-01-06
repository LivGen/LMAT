#!/bin/sh

if [ -z "$LMAT_DIR" ] ; then
   echo "Please set METAG_DIR environment variable to point to the directory where LMAT datafiles are stored"
   echo "Example syntax for bash/sh: export METAG_DIR=/data/here  csh: setenv METAG_DIR /data/here"
   exit 1
fi

## Should improve this
## Location of binaries
if hash rand_read_label >& /dev/null ; then
   bin_dir=
else
   bin_dir="$LMAT_DIR/../bin/"
   echo "Did not find rand_read_label in your path, assume LMAT binaries/scripts are here: $bin_dir"
fi

## platform specific support for hyperthreading
if [ -e /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus ] ; then 
   /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus
   export GOMP_CPU_AFFINITY=0-79
fi

## Assume the perm-je library is here
export LD_LIBRARY_PATH=$LMAT_DIR/../lib:$LD_LIBRARY_PATH


db_file=""
conv=$LMAT_DIR/m9.32To16.map
depthf="$LMAT_DIR/depth_for_ncbi_taxonomy.segment.pruned.dat"
taxtree=$LMAT_DIR/ncbi_taxonomy.segment.pruned.dat.nohl
rankinfo=$LMAT_DIR/ncbi_taxid_to_rank.txt
tax_histo_cnt=$LMAT_DIR/tcnt.m9.20.tax_histo

odir=.
read_len=0
read_range=0
num_bases=10000000000
min_sample_size=100
#create GC context specific null models, where GC content is evenly divided among this many bins (particuarly useful for handling larger eukaryotic genomes)
binsize=10
threads=80
tax_dir=""
usage="Generate random null model 

Usage: $0 options 

option list:
   --db_file=$db_file : taxonomy/kmer database 
   --read_len=$read_len : generate null model for a single read length
   --read_range=$read_range (ex: 50:1000:20) :  generate null models for a range of read lengths (e.g. from read legnth 50 to 10000 in intervals of 20)
   --num_bases=$num_bases (defualt): number of random reads to generate per thread
   --tax_histo_cnt=$tax_histo_cnt : observed k-mers associated with each taxid in the database
   --min_sample_size=$min_sample_size
   --odir=$odir : place output here
   --threads=$threads : number of threads to use
   --tax_dir=$LMAT_DIR

example usage:
$0 --db_file=$db_file --read_len=80
WARNING MUST MANUALLY SET SCRIPT FOR CORRECT OUTPUT FORMAT (OLD OR NEW) current=$SCRPT

"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
debug=0
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --debug*)
      debug=1;;
   --db_file=*)
      db_file=$optarg;;
   --read_len=*)
      read_len=$optarg;;
   --read_range=*)
      read_range=$optarg;;
   --num_bases=*)
      num_bases=$optarg;;
   --min_sample_size=*)
      min_sample_size=$optarg;;
   --tax_histo_cnt=*)
      tax_histo_cnt=$optarg;;
   --odir=*)
      odir=$optarg;;
   --threads=*)
      threads=$optarg;;
   --tax_dir=*)
      tax_dir=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

if [ -n "$tax_dir" ] ; then
   conv=$tax_dir/m9.32To16.map
   depthf="$tax_dir/depth_for_ncbi_taxonomy.segment.pruned.dat"
   taxtree=$tax_dir/ncbi_taxonomy.segment.dat.nohl
   rankinfo=$tax_dir/ncbi_taxid_to_rank.txt
   tax_histo_cnt=$tax_dir/tcnt.m9.20.tax_histo
fi

if [ $read_len = 0 ] &&  [ $read_range = 0 ] ; then
   echo "${usage}"
   exit 1
elif [ $read_len -gt 0 ]; then
   beg=$read_len
   end=$read_len
   int=100
else
   beg=`echo $read_range | awk '{ split($0,a,":"); print a[1] }'` 
   end=`echo $read_range | awk '{ split($0,a,":"); print a[2] }'` 
   int=`echo $read_range | awk '{ split($0,a,":"); print a[3] }'` 
   echo "what $read_range"
fi
#tot_reads=$[num_reads*threads]

db_file_name=`basename $db_file`
while [ $beg -le $end ] ; do
   read_len=$beg
   tmpval=$[num_bases/read_len]
   num_reads=$[tmpval/threads]
   echo "Create null model. Read_length=$read_len Reads per thread=$num_reads Total reads=$num_bases"
   dbname=`basename $db_file`
   oname=$dbname.$read_len.$num_bases.rl_output
   ofile=$odir/$oname
   logfile="$ofile.log"


   ${bin_dir}rand_read_label -w $rankinfo -f $conv -g $num_reads -i $read_len $vstr -e $depthf -p -t $threads -d $db_file -c $taxtree -o $ofile -f $conv   >& $logfile


   sfile=$ofile.rand_lst
   sname=`basename $sfile`
   if [ -e $sfile ] ; then
      ## try to roll up low read tax id calls into higher order taxids to avoid a lot of spurious strain calls
      merge_cnts.py $sfile $taxtree $rankinfo $min_sample_size $tax_histo_cnt $odir/null.bin.$binsize.$sname $binsize
      #cp $sfile $odir/null.bin.$binsize.$sname
      if [ $debug == 0 ] ; then
         gzip $odir/null.bin.$binsize.$sname
      fi
   else 
      echo "warning no $sfile found"
   fi
   let "beg += int"
done
ls -1 null.bin.*rand_lst.gz | perl -ne 'if(/\.(\d+)\.\d+\.rl_output/) { $t=$1-19; print "$t $_";}' | sort -k1n,1n > $db_file_name.null_lst.txt
echo "Manually copy *.rand_lst.gz and null_lst.txt to (LMAT_DIR) $LMAT_DIR"
#############################
#example:
#cp *rand_lst.gz $LMAT_DIR
#cp null_lst.txt $LMAT_DIR
#############################
