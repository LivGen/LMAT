# Outputs full fasta
def Readfile(fileName):
    #Create dictionary
    content=[]
    with open(fileName) as f:
        while True:
            #strip newline
            line=f.readline().rstrip("\n")
            if line.startswith(">"):
                id=line.replace(">seq.","") #number is left
                fragment=f.readline().rstrip("\n")
                #Add into directory
                content.append((id,fragment))
            elif not line:
                break
            else:
                #Exception for when reading the .fa file
                content.append(line)
    f.close()
    return content

def Assemble_genome(left,right,overlap):

    """Input is tuple with fragment id and fragment
        output is the string of both left-right united"""
    left_id=int(left[0])
    right_id=int(right[0])
    left=left[1]
    right=right[1]
    overlap=overlap

    #computes overlap between two fragments
    if left.endswith('N'):
        #This is made for the two strings to link together
        right='N'+right

    if right_id < left_id+3:
        for i in range(len(left)):
            if left[i:] == right[:len(left)-i]:
            #there is an overlap here: left[i:]
                return (right_id,left[0:i-1]+left[i:]+right[len(left)-i:])
    else:
    #if overlap doesnt exist add N's
        return (right_id,left+'N'*overlap)

#IMPORTS
import re
import sys
import os

####DESCRIPTION:####
#THIS SCRIPT CREATES THE FINAL FASTA PER SEQID


if __name__== "__main__":
    header_file=sys.argv[1] # input the header_file
    header_file=header_file.rstrip("\n")
    wd='/usr/mic/post1/metagenomics/ref_sets/fasta/01012015update/sinceMar2014'#depending the header file we do the wd
    folder=header_file.replace(".headers","")
    folder_dir=header_file.replace("headers","dir")#name of the directory to look out for the sequences
    with open(header_file) as hf:
        #change folder
        os.chdir(folder_dir)
        for header in hf:

        #get seqid using re
            match=re.search('\[sequence_id (\d+)\]',header)
            seqid=match.group(1)
            genome=""
            os.chdir(seqid)
            filename=seqid+"_final.fasta"
            if os.path.exists(filename)==False:
                g=open(filename,'a')
                v2_fasta=seqid+"_v2.fa.frag"
                if v2_fasta in os.listdir("."):
                    reads = Readfile(v2_fasta)
                    genome=reads[0]
                    for i in range(1,len(reads)):
                        genome=Assemble_genome(genome,reads[i],50)
                    genome=genome[1]
                else:
                    genome=Readfile(seqid+".fa")[0]

        #open a new txt file giant_fasta and append header and genome information in
        #g=open(os.path.join(wd,folder_dir,folder+"+.giant_fasta.txt"))
                g.write(header)
                g.write(str(genome))
                g.close()
            os.chdir('..')
    hf.close()
    
