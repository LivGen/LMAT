!# bin/bash

for i in $(cat run.ls); do
 cd $i;
/usr/gapps/kpath/lmat/LMAT-1.2.6/bin/rs.sh "$i".fa.frag;
 cd ../; 
done   
