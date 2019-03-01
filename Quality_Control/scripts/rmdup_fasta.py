#! /usr/local/bin/python

import pandas as pd
import sys

header_file=sys.argv[1] #2015_fasta.headers
fasta_file=sys.argv[2]
out='u_'+fasta_file

all=pd.read_table(header_file, sep="\n",header=None)
dups=all[all.duplicated()]
dups=dups[0].tolist()
reps=[0]*len(dups) #see if we already saw it
b=open(fasta_file,'r')
o=open(out,'a')

for line in b:
  line=line.rstrip() 
  if line in dups and reps[dups.index(line)] != 1:
      #find position of line in dups and increment reps
      #if rep is 1 means we already have it so we pass it.
      reps[dups.index(line)]+=1
      o.write(line+"\n")
 
  elif line in dups and reps[dups.index(line)] == 1:
      next(b) #skip its genome 
      	
  else:
      o.write(line+"\n")	 

b.close()
o.close()

