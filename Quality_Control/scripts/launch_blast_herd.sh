#! /bin/sh

#MSUB -l partition=herd
#MSUB -A genome
#MSUB -l nodes=1
#MSUB -l ttc=20
#MSUB -l walltime=72:00:00
#MSUB -q pbio-ng
#MSUB -N job.10002
#MSUB -V

wd=$(pwd) #the all_bacteria.dir for example and inside there must be a list of sequences
id_file=$wd/human_id
num_jobs=50
in_file=$wd/../virus.due #for original case search $wd/run.ls
while read i; do
  odir=$i/data.$i.fa.frag
  fna=$odir/$i.fa.frag.cleaned.db.humanid.9606.fna
  fasta=$fna.fasta

  #Extracts the 9606 reads -> output is the *.fna
 # if [! -e $fna]; then
 #	 /usr/gapps/kpath/lmat/LMAT-1.2.6/bin/pull_reads_mc.sh --idfile=$id_file --file_lst=$odir/tlst.allgenes.7-14.20.db.$i.fa.frag
 # fi

  #Cd-hit -> Output is the *.fasta

  if [! -e $fasta]; then
  	/usr/gapps/kpath/cd-hit/cd-hit -i $fna -o $fasta
  fi
  
  #Launches the blast
  /usr/bin/time -v /usr/gapps/kpath/ncbi-blast-2.2.27+/bin/blastn -db /p/lscratche/torres49/blastdb/db/nt -query $fasta -outfmt 7 -word_size 16 > $fasta.txt
  
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
done < $in_file 

