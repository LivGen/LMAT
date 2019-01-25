#! /usr/bin/python

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
   m=re.match("failed to find name entry for tid= (\d+)",line) 
   if m.group(1).isdigit()==True:    
      old_tid=m.group(1) 
      # get updated
      updated,otn=ncbi._translate_merged((old_tid,))
      new_tid=list(updated)[0] 

      print str(old_tid)+"\t"+str(new_tid)+"\n"
a.close()
 
