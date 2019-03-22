#! /bin/sh 
#MSUB -l nodes=2
#MSUB -l walltime=24:00:00
python rmdup_fasta.py 2015_fasta.headers 2015_full.fasta
