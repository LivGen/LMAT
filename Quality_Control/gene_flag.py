#This report counts how many genes are being identified per sequence and how many of those genes
#actually belong to the sequence being observed 


#imports
import pandas as pd
import os
import sys
import re

wd="/usr/mic/post1/metagenomics/ref_sets/fasta/01012015update/sinceMar2014/"
#input are the directorates ex. all_*.dir
db="allgenes.7-14.20.db"

for folder in sys.stdin:
    folder=folder.rstrip("\n")
    pwd=wd+folder
    match=re.search('all_(.*?).dir',folder)
    #i.e bacteria,virus,other,etc
    category=match.group(1)
    sum_file="summary_"+str(category)
    summary=pd.read_table(sum_file)
    summary["Sp_TaxID"]=[x.replace(",","") for x in summary["Sp_TaxID"]]
    results=[]
    #Change directory after reading
    os.chdir(folder)
   #the specie taxid is wrapped by commas

    for seqid in summary["SeqID"]:
        #get taxid,specieid and geneid
        ids=summary.loc[summary["SeqID"]==seqid,["TaxID","Sp_TaxID","Gen_TaxID"]]
        name=summary.loc[summary["SeqID"==seqid],"Name"]
        os.chdir(seqid)
        for gfile in os.listdir("data."+seqid+".fa.frag"):
            if gfile.endswith(".genesummary"):
                gene_summary=pd.read_table(gfile, header=None)
                genes=list(set(gene_summary[5]))
                spindex=[i for i,item in enumerate(gene_summary[3]) if item in ids]
                if len(spindex)==0:
                    count_sp=0
                else:
                    spgenes=list(set(gene_summary[5][spgenes]))
                    count_sp=len(spgenes)
                results.append([seqid,name,str(len(genes))],str(len(spgenes)),str(round(len(spgenes)/float(len(genes)),4))])

            else:
            #Append results to list: Seqid,name, number of genes counted,taxid correct genes, proportions
                results.append([seqid,name,"NA","NA","File Not Created"])

    results_df=DataFrame(results,columns=["seqid","name","count_genes","sp_genes","proportion"])
    results_df.to_csv("gene.results_"+category,sep=" ")
