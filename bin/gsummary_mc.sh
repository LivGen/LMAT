#!/bin/sh 
### default options
file_lst=""
min_kmers=30
prog=gsummary.pl
min_tax_score=0
min_score=0.10
num_threads=0
taxfile="$LMAT_DIR/gn_ref2.txt.gz"
usage="Generate a gene profile
Usage: $0 options

option list:
   --min_kmers=$min_kmers (default)
   --min_score=$min_score (default)
   --min_tax_score=$min_tax_score (default)
   --file_lst=$file_lst (default)
   --taxfile=$taxfile (default)
   --threads=$num_threads (default)
"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --min_kmers=*)
      min_kmers=$optarg;;
   --min_tax_score=*)
      min_tax_score=$optarg;;
   --min_score=*)
      min_score=$optarg;;
   --file_lst=*)
      file_lst=$optarg;;
   --taxfile=*)
      taxfile=$optarg;;
   --threads=*)
      num_threads=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

bindir=""
if ! hash $prog >& /dev/null; then
   #echo "Warning could not find $prog will try $LMAT_DIR/../bin"
   bindir="$LMAT_DIR/../bin/"
fi


jobstr=""
if [ $num_threads -gt 0 ] ; then
   jobstr="--jobs $num_threads"
fi
tfile="$file_lst.tmp.$$"
if [ -e $tfile ] ; then
   rm -f $tfile
fi
while read file ; do
   echo "$file $min_score $min_kmers $min_tax_score" >> $tfile
done < $file_lst

${pbin}parallel --load 150% --progress $jobstr -a $tfile ${bindir}$prog

lst=""
while read file ; do
   lst="$file.$min_score.$min_kmers.genesummary $lst"
   save="$file.$min_score.$min_kmers.genesummary"
done < $file_lst

out=`echo $save | sed 's/output[0-9]*\.out/output/'`
echo "$lst" | ${bindir}combine_gs.pl $taxfile | sort -k2nr,2nr > $out

lst=""
save=""
while read file ; do
   lst="$file.$min_score.$min_kmers.genesummary.min_tax_score.$min_tax_score $lst"
   save="$file.$min_score.$min_kmers.genesummary.min_tax_score.$min_tax_score"
done < $file_lst

out=`echo $save | sed 's/output[0-9]*\.out/output/'`
echo "$lst" | ${bindir}combine_gs.pl $taxfile | sort -k2nr,2nr > $out
