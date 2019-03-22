#!/bin/python
import sys
import re

input=sys.argv[1]
mapping_gi_file="taxonomy/nucl_wgs.accession2taxid"
gt_dir={}
at_dir={}


def Get_Accession(header):
	h=re.match(">([0-9A-Z]{0,12})\.\d",header)
	if h is not None:
		accession=h.group(1)
		if accession in at_dir.keys():
			taxid=at_dir[accession]
			return(taxid)	
	else:
		return(None) 

with open(mapping_gi_file,"r") as gi_f:
	for line in gi_f:
		info=line.rstrip().split()
		if info[0] not in at_dir.keys():
			at_dir[info[0]]=info[2] #This way we reduce the info stored
		gt_dir[info[3]]=info[2]
gi_f.close()

out=input+".multi"
out=open(out,"w")
with open(input,"r") as t:
	for inline in t:
		components=inline.rstrip().split()
		ti=components[0]
		gi=components[1]
		header=" ".join(components[3:])
		#If the taxid="0000" find through gi
		if ti=='0000' and gi != '-1':
			if gi in gt_dir.keys():
				components[0]=gt_dir[gi]
				out.write("\t".join(components[0:3])+"\t"+header+"\n") 
			else:
				print("Error no ti for gi:"+ str(gi)+"\n")
		elif ti=='0000' and gi=='-1': 
			#try withdraw information through accession
			taxid=Get_Accession(header)			
			if taxid != None:
				components[0]=taxid 
			out.write("\t".join(components[0:3])+"\t"+header+"\n")

		else:
			out.write(inline) #print same line 
t.close()
