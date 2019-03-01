#!/usr/local/bin/python

### IMPORTS
import ftplib
import os,sys
from time import sleep
import pandas as pd
import gzip

###Variables for job submission####
sleepTime=5
batchSize=100

###Method###
def getFile(ftp, filename):
    try:
        #Unzip file and then write to local directory
        ftp.retrbinary("RETR " + filename ,gzip.open(filename, 'wb').write)
    except:
        print "Error"

##### Input is list of accessions to download
##Please refer to confluence entrance:
##Example list of accessions can be found here:

file='ncbi_update.assembly'
 #second column is the accession number
downloadDir="/home/lebronaldea1/genbank_2017/ncbi_update"
file_df= pd.read_table(os.path.join(downloadDir,file),header=None)
file_ls=file_df[1].tolist()

#Open ftp
link='ftp.ncbi.nih.gov'
user='anonymous'
email='lebronaldea1@llnl.gov'
ftp=ftplib.FTP(link)
ftp.login(user,email)
#for each accession, divide the accession and create wd

for accession in file_ls:
    spl=accession.split("_")
    Dir='/genomes/all/'spl[0] #GCA folder or other
    projectID=spl[1].split(".")[0] #projectID.1
    #The folders are divided by splitting full project ID in 3 digits
    f1,f2,f3=projectID[:3],projectID[3:6],projectID[6:9]
    ### Search fasta files in directory and check creation dates
    fastaDir=os.path.join(Dir,f1,f2,f3)
    #There is a fourth folder which has the projectID and finishes in v1

    #Get list of all content and check which folder that starts with projectID is the last modified

    #get the fasta


    #get all modified times of fastas and download the newest one
    LastAccessTime=os.path.getmtime()

    #The fasta for download ends in
    fastaFile=f4+"_genomic.fna.gz"
    getFile(ftp,fastaFile)

    #how to deposit this fasta in the server?

ftp.quit()
