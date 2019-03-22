#!/bin/sh 

#################################################################################################################
### SCRIPT TO RUN THE content summary PIPELINE
### 
### The primary purpose of this procedure is to generate k-mer distribution plots for each species and plasmid
### to better evaluate relative abundance of individual organism calls. k-mer distribution plots can 
### be produced for any named rank, but the default is for species, plasmid and genus calls
###
### k-mer distribution plots are tradditionally used to estimate coverage for the single organism case
### and this becomes a fundamentally hard problem for the metagenomics case since
### the correct choice of k depends on two unknown variables
### relative abundance of the organism (presumably unknown) and the percentage of the organism's genome that is rank
### specific (this can potentially be calculated but is not currently done, and will likely continue to change
### as new organisms are sequenced).
###
### Therefore a range of k are computed for user evalutation. By default k-mer distributions are calculated for 
### k=8,10,12,14 and 17.
###
### A summary file is produced for the species level showing the top strain match (where available)
### and the "second peak" in the k-mer distribution plot, assumed to estimate coverage, or -1 when
### a second peak does not exist
###
### Script requires that the LMAT read label process was run 
###
### 2) Call content_summ to generate summary of organisms present using read_label output
###
#################################################################################################################
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
if hash content_summ >& /dev/null ; then
   bin_dir=

elif [ -f content_summ ] ; then
    bin_dir=./
elif [ `basename $PWD` == "nclass" ] ; then
    bin_dir="../apps/"
 else
   #echo "Could not find read_label in your path assume LMAT binaries/scripts are here: $bin_dir"
   bin_dir="$LMAT_DIR/../bin/"
fi
## Assume the perm-je library is here
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LMAT_DIR/../lib
#####################################################
## DEFAULT PARAMETER SETTINGS FOR CLASSIFICATION
#####################################################

## overwrite output
overwrite=0
## NCBI taxonomy tree given in LMAT format
taxtree="$LMAT_DIR/ncbi_taxonomy.segment.pruned.dat.nohl"
## human readable reporting of taxonomy lineage
taxfile="$LMAT_DIR/ncbi_taxonomy_rank.segment.pruned.txt"
## Used by content summarization to map taxids to rank (can this be merged with $species_map?)
rankval="$LMAT_DIR/ncbi_taxid_to_rank.pruned.txt"
# set to 1 for debugging only (too much output for large runs)
verbose=0
# specify directory to place output
odir=.
lst=""
fastsum_file=""
## plasmids that were not associated with a chromosome did not 
## get special taxids, this file identifies these plasmids
## so that during content summarization, the minimum coverage parameter for
## plasmids is correctly applied
xtra_plas_file="$LMAT_DIR/low_numid_plasmids.txt"

# run content summary on marker
usage="Run LMAT pipeline 
Usage: $0 options 

option list:
   --ilst=$lst: list of files containing LMAT's read label output 
   --filesum=$fastsum_file : LMAT's fastsummary output file
   --verbose=$verbose : Only used for debugging single read queries (too much output for larger datasets)
   --odir=$odir : Place output in this directory (defaults to current)
   --overwrite (default=$overwrite) : overwrite output file if it exists 
   --version : display version number and exit

example usage:
$0 --ilst=file_lst_of_lmat_output --filesum=*.fastsummary 

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
      lst=$optarg;;
   --filesum=*)
      fastsum_file=$optarg;;
   --odir=*)
      odir=$optarg;;
   --verbose)
      verbose=1;;
   --overwrite)
      overwrite=1;;
  --version)
      ${bin_dir}conten_summ -V ; exit;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done


vstr=""
if [ $verbose == 1 ] ; then
   vstr="-y"
fi

if [ ! $odir == '' ]; then
    odir="$odir/"
fi


if [ -e $fastsum_file ] ; then
   fastsum_ofile_name=`basename $fastsum_file`
   fastsum_ofile="${odir}$fastsum_ofile_name"
   sumofile="$fastsum_ofile.summ"
   logfile="$sumofile.log"
   kcnt="8,10,12,14,17"
   hstr=""
   if [ ! -e $sumofile ] || [ $overwrite == 1 ] ; then
      rank_check="plasmid,species,genus"
      echo "Summary process $fastsum_file [overwrite=$overwrite (1=yes, 0=no)] [outputfile=$sumofile] [k-values=$kcnt]"
      /usr/bin/time -v ${bin_dir}content_summ $hstr -p $xtra_plas_file -c $taxtree -l $fastsum_file -k $kcnt -f $lst -r $rankval -a $rank_check -o $sumofile >& $logfile
      ofile="$fastsum_ofile.ordered"
      ${bin_dir}summary.py $sumofile $rankval $fastsum_file $xtra_plas_file $ofile $rank_check
   else 
      echo "$sumofile exists, use --overwrite to overwrite file"
   fi
else 
   echo "Error requires LMAT's *fastsummary file as input"
fi
