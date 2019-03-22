#! /usr/local/bin/python
#Dayanara Lenron
#01/19
#This script removes duplicate entries in a fasta file.

import pandas as pd
import sys

fasta_file=sys.argv[1] #2015.fasta

out='nodup_'+fasta_file

all_headers=[]
b=open(fasta_file,'r')
o=open(out,'a')

for line in b:
  if len(line.strip())==0:
	continue
  line=line.rstrip()
#Specify if it is duplicated
  if line.startswith(">") and (line not in all_headers):	
	all_headers=all_headers.append(line)
	dup=0
  elif line.startswith(">") and line in all_headers:
	dup=1
#Write lines if not duplicated
  if dup==0:
      o.write(line+"\n")
  else:
      next(b) #skip its genome 

b.close()
o.close()

