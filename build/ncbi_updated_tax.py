#! /usr/bin/python
#Dayanara Lebron 
#This python script gets updated taxids fro

from sys import *
import os
import re
import ete3
from ete3 import NCBITaxa

ncbi=NCBITaxa()

a=open(argv[1]) # this is the file that has all the errors
a.readline()
for line in a:
   line=line.rstrip()
   if line.startswith("failed"):
	m=re.match("failed to find name entry for tid= (\d+)",line) 
	old_tid=m.group(0)
   elif line.startswith("weird:"):
	comp=line.split()
	old_tid=comp[2]
   else:
	continue
	
   if old_tid.isdigit()==True:     
      # get updated
      updated,otn=ncbi._translate_merged((old_tid,))
      new_tid=list(updated)[0] 

      print str(old_tid)+"\t"+str(new_tid)+"\n"
a.close()
 
