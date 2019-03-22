#!/bin/sh

#####################################
### SCRIPT TO RUN THE Gene labeling PIPELINE
###
### Steps:
###
### Call gene_label to assign genes to each read and count the reads assigned to each gene 
###
#####################################
if [ -z "$LMAT_DIR" ] ; then
   echo "Please set LMAT_DIR environment variable to point to the directory where LMAT datafiles are stored"
   exit 1
fi

## Some environments require explicit enabling of hyperthreading
## Other environments may already enable this
if [ -e /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus ] ; then
   /collab/usr/global/tools/mpi/utils/hyperthreading/enable_cpus 
   export GOMP_CPU_AFFINITY=0-79
fi

## Should improve this 
## Location of binaries
if hash gene_label >& /dev/null ; then
   bin_dir=
elif [ -f gene_label ] ; then
    bin_dir=./
elif [ `basename $PWD` == "nclass" ] ; then
    bin_dir="../apps/"
 else
   bin_dir="$LMAT_DIR/../bin/"
fi
## Assume the perm-je library is here
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LMAT_DIR/../lib




#####################################
## default parameter settings
#####################################
fs_file=""
overwrite=0
min_tax_score=0
## Memory mapped gene database file
genedbfile=""
## Used by gene_label to assign human readable names to genes
genefile="$LMAT_DIR/gn_ref2.txt.gz"
## ignore reads with less valid k-mers than this value
min_read_kmer=30
rank_report="species,plasmid,genus"
## minimum percentage of k-mer matches from a read found in a gene before a call can be made
gene_score=0.1

## Additional user input default settings
# set to 1 for debugging only (too much output for large runs)
verbose=0
# specify directory to place output
odir=.
## minimum number of valid k-mers  present in read needed before considering the read for gene labeling
num_gene_kmers=20
lstr=""

# in some cases you might not want to bother summarizing the human reads
# turn on if there are a huge # of human reads and you want to save time

# run content summary on marker
usage="Run LMAT pipeline 
Usage: $0 options 

option list:
   --ilst=$lstr : File list LMAT read_label output
   --db_file=$genedbfile : Memory mapped gene database file
   --verbose=$verbose : Only used for debugging single read queries (too much output for larger datasets)
   --odir=$odir : Place output in this directory (defaults to current)
   --filesum=$fs_file : fastsummary file generated from read_label. When file is specified (optional)
                        a post processing script is automatically invoked to provide an updated taxonomy report, which includes additional information on gene content
                        reports percentage of rRNA for each taxon and number distinct gene ids detected and number of reads mapped to these genes
                        this feature is meant to give additional information on how much of the genome was used for taxnomic identification
   --overwrite (default=$overwrite) : overwrite output file if it exists 
   --min_read_kmer (default=$num_gene_kmers) : minimum number of valid k-mers present in read needed for analysis
   --min_gene_score (default=$gene_score) : minimum percentage of k-mers matching to reference gene (for gene summary step)
   --min_tax_score (default=$min_tax_score) : minimum score for matched tax id (used for tracking rRNA assigned to specific tax ids)
   --rank_report=$rank_report : read binning for different ranks (user provides a comma separated list of ranks). plasmid is treated as a separate rank
   --version : report version number and exit

example usage:
$0 --db_file=$genedbfile --ilst=run_lmat_output_file_lst.lst 

"

if test $# = 0; then
   echo "${usage}"
   exit 1
fi

while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --ilst=*)
      lstr=$optarg;;
   --odir=*)
      odir=$optarg;;
   --min_read_kmer=*)
      num_gene_kmers=$optarg;;
   --db_file=*)
      genedbfile=$optarg;;
   --filesum=*)
      fs_file=$optarg;;
   --verbose)
      verbose=1;;
   --overwrite)
      overwrite=1;;
   --min_gene_score*)
      gene_score=$optarg;;
   --rank_report)
      rank_report=$optarg;;
   --min_tax_score)
      min_tax_score=$optarg;;
  --version)
      ${bin_dir}read_label -V ; exit;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done

min_gene_read=0

vstr=""
if [ $verbose == 1 ] ; then
   vstr="-y"
fi

if [ ! $odir == '' ]; then
    odir="$odir/"
fi

## assign gene names
if [ $genedbfile ] ; then
      genedbname=`basename $genedbfile`
      query_file_name=`basename $lstr`
      qstr="-l $lstr"
      vstr=""
      ## note need to fix verbose setting to get here (if needed)
      if [ $verbose == 1 ] ; then
         vstr="-y"
      fi
      genofile="${odir}$query_file_name.$genedbname.gl_output"
      logfile="${odir}$query_file_name.$genedbname.gl_output.log"
      res=$genofile.$gene_score.$num_gene_kmers.genesummary
      res2=$genofile.$gene_score.$num_gene_kmers.genesummary.min_tax_score.$min_tax_score
      if [ ! -e $res ] || [ $overwrite == 1 ] ; then
         ${bin_dir}gene_label $vstr -b $min_tax_score -q $num_gene_kmers -x $gene_score -p $qstr -d $genedbfile -o $genofile -g $genefile >& $logfile
         cat $res | sort -k1gr,1gr > tmp.$$
         mv tmp.$$ $res 
         min_map_reads=10
         if [ ! -z $fs_file ] ; then
            ${bin_dir}fsreport.py $fs_file $rank_report $odir $res2 $min_map_reads
         fi
      fi
else
   echo "Could not localte database file: $genedbfile"
   exit 1
fi
