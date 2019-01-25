#!/bin/sh 
### default options
file_lst=""
num_jobs=200
usage="Manage submission of over 200 jobs to adhere to 200 job submission limit
Usage: $0 options

option list:
   --num_jobs=$num_jobs (default)
   --file_lst=$file_lst (default)
"
if test $# = 0; then
   echo "${usage}"
   exit 1
fi
while test -n "${1}"; do
   opt=${1}
   optarg=`expr "x$opt" : 'x[^=]*=\(.*\)'`

   case $opt in
   --num_jobs=*)
      num_jobs=$optarg;;
   --file_lst=*)
      file_lst=$optarg;;
   *)
      echo "Unrecognized argument [$opt]"
      echo "${usage}"
      exit 1
   esac
   shift
done
user=`whoami`
while read cmd ; do
   #sbatch -o log.$file --time=24:00:00 pcmd.sh $file
   #sbatch -o log.$file --time=24:00:00 pcmd.sh $file
   $cmd
   doWait=1
   while [ $doWait -eq 1 ] ; do
      runcnt=`squeue -u bioinf | grep bioinf | wc -l`
      if [ "$runcnt" -ge "$num_jobs" ] ; then
         doWait=1
         sleep 5
      else
         doWait=0
      fi
   done
done < $file_lst
