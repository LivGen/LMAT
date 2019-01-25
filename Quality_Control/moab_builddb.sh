#!/bin/tcsh
#MSUB -l walltime=48:00:00
#MSUB -o status.dbbuild
#MSUB -l nodes=2
#MSUB -q pbio-ng 

#size=900
#wdir=.
#cutoff=200
#ramdisk=/l/ssd
#oname='18102017'

#devdir=/usr/gapps/kpath/lmat/utils/dbbuild/

srun /usr/gapps/kpath/lmat/utils/dbbuild/build2_db.sh 18102017update 900 200 . /l/ssd
echo 'Done'

