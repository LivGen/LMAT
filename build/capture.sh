#!/bin/sh


odir=$1

while :
do 
   pid=`ps -C read_label | grep read_label | cut -d" " -f1`
   if [ "$pid" ]; then
      cat /proc/meminfo > $odir/save.meminfo.$pid
      cat /proc/$pid/status > $odir/save.status.$pid
   fi
   sleep 120
done
