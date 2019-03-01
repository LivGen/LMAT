#!/usr/bin/env python

import pandas as pd
import re
import sys
import ete3
from ete3 import NCBITaxa

ncbi=NCBITaxa()

# Read reference file
rfl='assembly_summary_genbank.txt' #prior fix by removing first line and #'s'
reference=pd.read_table(rfl)

for header in sys.stdin: #change
    header=header.rstrip()
    #Retrieve taxid information
    mtch=re.match("(.*?)\[tax\_node\_id (\d+)\]", header)
    try:
        tid=mtch.group(2)
    except:
        #If there is no match
        print("error: no taxid for ",header,"\n")
        continue

    #Get info from reference about taxid
    try:
        found=reference[reference["taxid"]== int(tid)]
	tid=found["species_taxid"]

    except:
        #this might be an old taxid or a merged taxid; get new taxid
        updated, otn=ncbi._translate_merged((tid,))
       	tid=updated
        found=reference[reference["taxid"]== int(tid)]
	tid=found["species_taxid"]

        #Maybe make exception for plasmids with taxid 1000
    obs=found.shape[0]

    #Get accession version
    #if 'strain' in header:
    #    try:
    #	#try to match accession for strain information and match
    #       t=re.match("(.*?) strain (.*?),",header)
    #       stra='strain='+t.group(2)
    #	except:
    #	   print("error in header getting strain info:",header,"\n")
	   	
    if int(obs) == 1:
        accession=found["assembly_accession"]
      #Multiple strains 
    if len(accession)>2:
	#Get info for specific strain
        accession=accession[1]
	#That strain has multiple release dates/ duplicates
        #if specific_strain.shape[0] > 1:
            #just get the first line
            #it can be done per seq_rel_date in this manner but some dont match even then
            #Example: tid=1000570; release: 2011 and 2013 sequence from 2015.
         #   accession=specific_strain.iloc[0]["assembly_accession"]
        #else:
         #    accession=specific_strain["assembly_accession"]
    print(accession+"\n")
    print(tid+"\n") 	 
    print(accession+"\t"+tid+"\n")

