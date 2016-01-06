#!/bin/sh 
### default options
file_lst=""
odir=.
taxfile=$LMAT_DIR/ncbi_taxonomy_rank.segment.pruned.txt
num_threads=0 ## default to number of cores on machine
sstr=""
idfile=""
usage="Build a taxonomy identifier list based on standard natural language taxonomy nomenclature
Usage: $0 options

option list:
   --taxfile=$taxfile : include reads with >= this number of valid k-mers - useful when N masking is used
   --file_lst=$file_lst (default)  : list of *.fastsummary files with LMAT output. Assumption is there are *.fastsummary.krona files (default output format for LMAT)
   --idfile=$idfile (default) : place output in this file
   --sstr=$sstr  : search string,  examples "kingdom,Virus"  or "genus,Bacillus", see $taxfile for format 
Example
    build_taxid_lst.sh --file_lst=fs.flst --sstr=superkingdom,Viruses --idfile=none

"

if test $# = 0; then
   echo "${usage}"
   exit 1
fi

bindir=""
pbin=""
if ! hash parallel >& /dev/null; then
   pbin="$LMAT_DIR/../bin/"
fi

if ! hash build_taxid_lst.pl >& /dev/null; then
   bindir="$LMAT_DIR/../bin/"
fi
   
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --taxfile=*)
      taxfile=$optarg;;
   --file_lst=*)
      file_lst=$optarg;;
   --idfile=*)
      idfile=$optarg;;
   --sstr=*)
      sstr=$optarg;;
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
   ofile=$file.idtxt
   echo "${bindir}build_taxid_lst.pl $taxfile $file $ofile $sstr | ${pbin}parallel " >> $tfile
done < $file_lst
${pbin}parallel --load 150% --progress $jobstr -a $tfile 

rm -f $tfile
