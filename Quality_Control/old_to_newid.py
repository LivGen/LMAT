#!/bin/python 
#This script takes as input the oldids gotten from the remap and gets the new ids 
import ete3
import pandas as pd
from ete3 import NCBITaxa
ncbi = NCBITaxa()

import sys

f=open(sys.argv[1])
ids=[i.replace(")","").replace("'","").rstrip().strip() for i in f]

new_id=[ncbi.get_lineage(int(id))[-1] for id in ids] # a message will show that the id has been merged into another

#Could be done with pandas but some LC do not have it installed
o=open("new_old","w")
outstring=zip(ids,new_id)
for i in outstring:
	o.write(" ".join(x for x in line)+"\n")
o.close()
f.close()
