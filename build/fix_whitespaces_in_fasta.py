#!/bin/python
#Dayanara Lebron Aldea
#01/19
#This script removes whitespaces in the input fasta file such that a genome 
#appears in one line instead of fragmented.

import sys
fasta_file=sys.argv[1]

f=open(fasta_file)
o=open("no_trail."+fasta_file,"w")
genome=""
for line in f:
#	print(line)
	if not line.startswith(">"):
		#Take the whitespace out of the fragment 
		line=line.rstrip()
		genome+=line #Concatenate genome fragments
	else:	
		if genome !="":
			o.write(genome+"\n") #Write out concatenated genome
		o.write(line) #This is a header
		genome="" #Reset for concatenation
o.write(genome+"\n") #So that the last genome can be written
o.close()
f.close()

			
