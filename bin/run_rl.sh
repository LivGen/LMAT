#!/bin/sh 

#####################################
### SCRIPT TO RUN THE LMAT PIPELINE
###
### Steps:
###
### 1) Call read_label to taxonomically label all reads, and count reads assigned to each label
### 2) call tolineage.py generates human readable taxonomy lineage, can be input to Krona, which is then run if the binary is found, to produce an html file
### 3)  Call fsreport to report at genus, species and strain levels, plus a plasmids report
### 4) Generate html output from genus and species/strain report
###
#####################################
if [ -z "$LMAT_DIR" ] ; then
   echo "Please set LMAT_DIR environment variable to point to the directory where LMAT datafiles are stored"
   exit 1
fi

## Some environments require explicit enabling of hyperthreading
## Other environments may already enable this
#if [ -e /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus ] ; then
#   /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus 
#   export GOMP_CPU_AFFINITY=0-79
#fi

## Should improve this 
## Location of binaries
if hash read_label >& /dev/null ; then
   bin_dir=

elif [ -f read_label ] ; then
    bin_dir=./
elif [ `basename $PWD` == "nclass" ] ; then
    bin_dir="../apps/"
 else
   bin_dir="$LMAT_DIR/../bin/"
fi
## Assume the perm-je library is here
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LMAT_DIR/../lib
overwrite=0
#####################################################
## DEFAULT PARAMETER SETTINGS FOR CLASSIFICATION
#####################################################
## Memory mapped taxonomy database file
dbfile=""

## NCBI taxonomy tree given in LMAT format
## this version removes most of the intermediate human lineage nodes
## to save compute time
taxtree="$LMAT_DIR/ncbi_taxonomy.segment.pruned.dat.nohl"
## content_caller_qsum.py  uses the non pruned version
noprune_taxtree=$LMAT_DIR/ncbi_taxonomy.segment.pruned.dat
## Stores tree depth of each taxnomy node
depthf="$LMAT_DIR/depth_for_ncbi_taxonomy.segment.pruned.dat"
## identify the rank for each taxid - uses a keyword "strain" to identify ranks below species
rankval="$LMAT_DIR/ncbi_taxid_to_rank.pruned.txt"
## Stores human readable lineage information on each taxonomy node
taxfile="$LMAT_DIR/ncbi_taxonomy_rank.segment.pruned.txt"
## For content summarization, and read counting, ignore reads with scores below this threshold
min_score=0

## Should be deprecated with version 1.2.3 (it increases the human tax scores by 1 standard deviation)
hbias=0

## The higher the number the more conservative the read label call
## This value specifies how much higher (in standard deviation units) the score of the assigned label must be
## relative to the competing taxonomic calls
sdiff=1.0  

## ignore reads with less valid k-mers than this value
min_read_kmer=30

# number of threads
threads=80
## unless otherwise specified will default to using
## all available cores
if [ -e /proc/cpuinfo ] ; then
   threads=`grep -c processor /proc/cpuinfo`
fi
# set to 1 for debugging only (too much output for large runs)
verbose=0
# Must be specified by user on command line
query_file=
# specify directory to place output
odir=.  

# run content summary on marker
usage="Run LMAT pipeline 
Usage: $0 options 

option list:
   --db_file=$dbfile : Memory mapped database file
   --query_file=$query_file : Metagenomic reads in fasta format
   --threads=$threads : Number of threads to use (will default to using all available cores) 
   --nullm=$nullm : File containing the list of null models
   --verbose=$verbose : Only used for debugging single read queries (too much output for larger datasets)
   --sdiff=$sdiff : Scoring differential. Tax ids with scores at or above the maximum score - sdiff*STDEV are considered
   --hbias=$hbias : For human samples where human DNA concentration is high, human taxid score + hbias*STDEV, hbias > 0 may be conveniant
                    (with latest database this parameter will likely no longer be needed) 
   --odir=$odir : Place output in this directory (defaults to current)
   --min_score=$min_score : minimum score assigned to read for it to be included in binning
   --overwrite (default=$overwrite) : overwrite output file if it exists 
   --min_read_kmer (default=$min_read_kmer) : minimum number of valid k-mers present in read needed for analysis
   --prune_thresh=X : threshold of maximum taxonomy IDs allowed per k-mer (default is to use default settings of prepruned database)
   --version : print version number and exit

example usage:
$0 --db_file=$dbfile --query_file=query.fna 

"

if test $# = 0; then
   echo "${usage}"
   exit 1
fi

while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`
   case $opt in
   --min_score=*)
      min_score=$optarg;;
   --hbias=*)
      hbias=$optarg;;
   --odir=*)
      odir=$optarg;;
   --db_file=*)
      dbfile=$optarg;;
   --min_read_kmer=*)
      min_read_kmer=$optarg;;
   --query_file=*)
      query_file=$optarg;;
   --prune_thresh=*)
      PTHRESH=$optarg;;
   --threads=*)
      threads=$optarg;;
   --nullm=*)
      nullm=$optarg;;
   --sdiff=*)
      sdiff=$optarg;;
   --verbose)
      verbose=1;;
   --overwrite)
      overwrite=1;;
   --version)
      ${bin_dir}read_label -V ; exit;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

fastq="false"

if [ ! -e $query_file ] ; then
   echo "Error $query_file not found"
   exit 0
fi

first=`head -n1 $query_file | awk '{print substr($1,0,1)}'`
fqchar='@'
    
if [ $first == $fqchar ] ; then
    fastq="true"
else 
    fachar='>'
    if [ $first != $fachar ] ; then
	echo "Error: $query_file must be in fasta or fastq format"
	exit 1
    fi
fi

if [ $fastq == "true" ] ; then 
    fastqstr="-q"
fi


if [ ! -e $db_file ] ; then
   echo "Error need to supply a markery library or full database file"
   exit 0
fi

query_file_name=`basename $query_file`

vstr=""
if [ $verbose == 1 ] ; then
   vstr="-y"
fi

if [ ! $odir == '' ]; then
    odir="$odir/"
fi

db=$dbfile
dbname=`basename $db`
if test -z $PTHRESH; then
  rlofile="${odir}$query_file_name.$dbname.lo.rl_output" 
else
  rlofile="${odir}$query_file_name.$dbname.$PTHRESH.lo.rl_output"	 
fi

logfile="$rlofile.log" 
tidmap="$LMAT_DIR/m9.32To16.map"
## File giving a list of null models - assumes this specific naming convention

fstr="-f $tidmap"
if [ ! -z $PTHRESH ]; then
   pstr="-g $PTHRESH -m $LMAT_DIR/numeric_ranks"
fi
rprog=${bin_dir}read_label
use_min_score=$min_score

if [ -z $nullm ] ; then
   echo "Using default null model list file for $dbname"
   nullm=$LMAT_DIR/$dbname.null_lst.txt
   nullmstr="-n $nullm"
elif [ "$nullm" == "no" ] ; then
	  echo Not using null model files
	  nullmstr=""
elif [ -f $nullm ] ; then
	  nullmstr="-n $nullm"
else
	  echo "Please provide valid null models list file"
	  echo "no file found at $nullm"
	  exit 
fi

fastsum_file="$rlofile.$use_min_score.$min_read_kmer.fastsummary"
echo "Create output file from input parameters: $fastsum_file"
if [ ! -e $db ] ; then
   echo "Could not find db: $db"
   exit 0
fi
if [ ! -e $fastsum_file ] || [ $overwrite == 1 ] ; then
   ostr=""
   if [ $overwrite == 1 ] ; then
      ostr="(overwrite existing files)"
   fi
   echo "Process $query_file outputfile=$fastsum_file $ostr"

   /usr/bin/time -v $rprog $fstr $pstr -u $taxfile -w $rankval -x $use_min_score -j $min_read_kmer -l $hbias -b $sdiff $vstr $nullmstr -e $depthf -p -t $threads -i $query_file -d $db -c $taxtree -o $rlofile $fastqstr >& $logfile
   min_reads=1
   if [ ! -e $fastsum_file ] ; then
      echo "Error, did not create a fastsummary file [$fastsum_file]"
      exit 0
   fi
   min_num_reads=10 ## minimum abundance cutoffs to help with Krona plots
   min_avg=0 ## minimum average read score (off)
   ${bin_dir}tolineage.py $taxfile $fastsum_file $fastsum_file.lineage $min_num_reads $min_avg
   ${bin_dir}fsreport.py $fastsum_file plasmid,species,genus $odir


#   python ${bin_dir}genusspecies2html.py ${odir}$fastsum_file.species ${odir}$fastsum_file.genus $taxfile > ${odir}$fastsum_file.html



   if hash ktImportText > /dev/null 2>&1 ; then
      ktImportText $fastsum_file.lineage -o $fastsum_file.lineage.html
   fi
else 
   echo "Warning, $fastsum_file exists, set --overwrite to overwrite existing file"
fi
