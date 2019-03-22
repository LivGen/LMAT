#!/bin/sh 
### default options
file_lst=""
min_kmers=15
odir=.
min_score=1
num_threads=0 ## default to number of cores on machine
idfile=""
prog=pull_reads.pl
usage="Generate bam file
Usage: $0 options

option list:
   --kmers=$min_kmers : include reads with >= this number of valid k-mers - useful when N masking is used
   --min_score=$min_score (default) : minimum score to use
   --file_lst=$file_lst (default)  : list of files with LMAT output. Multiple files are assumed to run search in parallel
   --idfile=$idfile (default) : list of taxids to pull out. Each line is a space separated list of taxids. One file of reads is created for each line of taxids specified in the file
                              : For example a file containing
                              : 1236 1149864
                              : 1030145
                              : LowScore 0
                              : NoDbHits
                              : will create 4 files one containing reads assigned taxids with value 1236 or 1149864, the second file will contain reads assigned 1030145, another file containing reads with scores less than 0 and the last file contains reads with no match to the database
                              : The first id is included in the filename to differentiate files so it is currently assumed to be unique.
   --threads=$num_threads (default) : default is to use all cores - you can customize the call to gnu parallel to limit resource use.
   --odir=$odir (default) : output directory

"

if test $# = 0; then
   echo "${usage}"
   exit 1
fi

## Some environments require explicit enabling of hyperthreading
## Other environments may already enable this
if [ -e /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus ] ; then
   /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus
   export GOMP_CPU_AFFINITY=0-79
fi

bindir=""
if ! hash $prog >& /dev/null; then
   echo "Warning could not find $prog will try $LMAT_DIR/../bin"
   bindir="$LMAT_DIR/../bin/"
fi

pbin=""
if ! hash parallel >& /dev/null; then
   echo "Warning could not find GNU parallel will try $LMAT_DIR/../bin"
   pbin="$LMAT_DIR/../bin/"
fi
   
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --kmers=*)
      min_kmers=$optarg;;
   --min_score=*)
      min_score=$optarg;;
   --file_lst=*)
      file_lst=$optarg;;
   --idfile=*)
      idfile=$optarg;;
   --odir=*)
      odir=$optarg;;
   --threads=*)
      num_threads=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

jobstr=""
if [ $num_threads -gt 0 ] ; then
   jobstr="--jobs $num_threads"
fi

tfile="$file_lst.tmp.$$" 
if [ -e $tfile ] ; then
   rm -f $tfile
fi
while read file ; do
   echo "$file $idfile $min_score $min_kmers $odir" >> $tfile
done < $file_lst

${pbin}parallel --load 150% --progress $jobstr -a $tfile ${bindir}$prog 
while read idlst ; do
   taxid=`echo $idlst | awk '{ print $1 }'`
   first=0
   while read file ; do
      bname=`basename $file`
      idname=`basename $idfile`
      ofile="$odir/$bname.$idname.pulled.$taxid"
      if [ -e $ofile ] ; then
         ## take this opportunity to shorten the filename
         mergefile=`echo $ofile | sed -r 's/lo.rl_output[0-9]+.out.//'`
         mergefile=`echo $mergefile | sed -r 's/pulled.//'`
         mergefile="$mergefile.fna"
         if [ $first -eq 0 ] ; then
            cat $ofile > $mergefile
         else 
            cat $ofile >> $mergefile
         fi
         first=1
         rm -f $ofile
      else 
         echo "Warning could not find $ofile"
      fi
   done < $file_lst
done < $idfile

rm -f $tfile
