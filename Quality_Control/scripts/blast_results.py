############################
#Argument #1: List of blast files:
############################

# Module Imports
import pandas as pd
import re
import sys
import numpy as np
import collections
import os

#READ ALL OF THESE FILES ONCE
wd="/usr/mic/post1/metagenomics/ref_sets/fasta/01012015update/sinceMar2014/"
wd_tax=wd+"../taxonomy/2017tax/"
#Human Genes accession
#to check if gene is human

human_genes=pd.read_table(wd_tax + "human_genes_gb",header=None)

#accession to taxid
acc_to_tax=pd.read_table(wd_tax + "nucl_gb.accession2taxid")

#domain info
#to check if taxid is from an eukaryota
taxid_to_k=pd.read_table(wd_tax + 'categories.dmp', header=None, names=["Domain","sptaxid","taxid"])

#Get taxonomy information
taxnames=pd.read_table(wd_tax+'names.dmp', header=None, names=["taxid","","name","","variant","","name_class",""])
taxnames=taxnames[["taxid","name","name_class"]].loc[taxnames["name_class"]=="scientific name"]

#Variables
threshold= 1e-10

#Getting location of fasta files
pwd = os.getcwd()
if 'Bact' in pwd:
    fasta_loc= wd+"all_bacteria.dir"
    catalog_loc=wd+"Data/catalog.bacteria"
elif 'Virus' in pwd:
    fasta_loc= wd+"all_virus.dir"
    catalog_loc=wd+"Data/catalog.virus"
elif 'Fungi' in pwd:
    fasta_loc= wd+"all_fungi.dir"
    catalog_loc=wd+"Data/catalog.fungi"
elif 'Other' in pwd:
    fasta_loc= wd+"all_other.dir"
    catalog_loc=wd+"Data/catalog.other"
else:
     fasta_loc= wd+"all_protozoa.dir"
     catalog_loc=wd+"Data/catalog.protozoa"

#Create Human_Report for amount of sequences identified by BLAST
human_report=[]

#Start reading the incoming blast files
for file in sys.stdin:
    file=file.rstrip("\n")
    seqid=file.split(".")[0]
    print(seqid)
    catalog=pd.read_table(catalog_loc, header=None, names=["name","seqid","taxid"])
    catalog["seqid"]=catalog["seqid"].astype(str)
    seq_taxid=catalog.loc[int(np.where(catalog["seqid"]==seqid)[0]),"taxid"]
    seq_name=catalog.loc[int(np.where(catalog["seqid"]==seqid)[0]),"name"]
    info=[]
    frag=0 #Number of fragments
    with open(file) as fp:
        for line in fp:
            match=re.search(r'^seq.(\d+)',line)
            if match:
                frag+=1
                seq_frag=match.group(1)
                components=line.split('\t')
                match2=re.search('gi\|\d+\|(.*?)\|(.*?)\|',components[1])
                accgene=match2.group(2)
                evalue=components[10]
                info.append([int(seq_frag),accgene,float(evalue)])

    df=pd.DataFrame(info,columns=['frag','accession','evalue'])
                #Keep significant e-values and examine
    df=df[df['evalue'] <= threshold]

    #check if any true then there are some human genes
    check=df['accession'].isin(human_genes[0])
    if(any(check)==True):
        delete=list(set(df["frag"][check]))
    else:
        delete=[]

    #Get  the taxid
    set_acc=set(df["accession"])
    indices=[i for i, item in enumerate(acc_to_tax['accession.version']) if item in set_acc]
    gene_taxid=list(acc_to_tax['taxid'][indices].astype(int))
    acc_v=list(acc_to_tax['accession.version'][indices])
    # Get domain for each taxid
    index_dict = dict((value, idx) for idx,value in enumerate(taxid_to_k['taxid']))

    indices_dom=[index_dict.get(x) for x in gene_taxid]

    #if index_dict.get(x) # is not None
    gene_domain=taxid_to_k['Domain'][indices_dom]
    gene_info=pd.DataFrame(np.column_stack((acc_v, gene_taxid,gene_domain)),columns=['gene.acc','taxid','Domain'])

    #If the seq_taxid is one of the matches then we can erase from the list, we are trying to report the false positives
    check=seq_taxid in gene_info["taxid"]
    if check:
        gene_info=gene_info.loc[gene_info["taxid"] != seq_taxid]

    #Remove the Eukaryota Fragments from the original file
    euk_frag=list(set(df[df["accession"].isin(gene_info["gene.acc"].tolist())]["frag"].tolist()))

    rmv_total=delete+euk_frag #add the human fragments

    #how many fragments were labeled as that gene
    pairs=set(list(zip(df["frag"], df["accession"])))
    counts=dict()
    for pair in pairs:
    #pair[1]=gene_name
        if pair[1] not in counts:
            counts[pair[1]]=1
        else:
            counts[pair[1]]+=1

    gene_info["counts"]=[counts[gene] for gene in gene_info["gene.acc"]]

    # get taxname for taxid of genes.
    gene_info["tax_name"]=[taxnames["name"][taxnames["taxid"]==taxid].tolist()[0] for taxid in gene_info["taxid"]]

    #Print report
    #Report file - ?Report removed sequences?
    gene_info.to_csv(seqid+"_blast.report.txt", index=None, sep=' ', mode='a')

    #Create fasta_v2
    #Re-write the frag file without the eukaryota and human fragments
    # Write human_report on all seqid

    if frag==0:
	    prop=0
	    prop_euk=0
    else:
	    prop=str(round(len(delete)/float(frag),4))
   	    prop_euk=str(round(len(euk_frag)/float(frag),4))

    #Reporting: Seq_id,number of fragments and proportion of overall fragments
    human_report.append([seqid,seq_name,frag,len(rmv_total),prop, prop_euk])

    #If there are fragments to remove then create the -v2 fasta
    if len(rmv_total)> 0:

        f=open(fasta_loc+"/"+seqid+"/"+seqid+'.fa.frag',"r")
        p=open(fasta_loc+"/"+ seqid+"/"+seqid+'_v2.fa.frag',"w")

        for line in f:
            match=re.search('^>seq.(\d+)',line)
            if match and int(match.group(1)) in rmv_total:
                next(f)
        #next #skips sequence line
            else:
                p.write(line)
    	f.close()
    	p.close()
df=pd.DataFrame(human_report, columns=["seqid","seq_name","frag_in_fasta","removed-frag",
"prop-Human","prop-Euk"])
df.to_csv("Human_Report.txt",index=None,sep=" ", mode='a')
