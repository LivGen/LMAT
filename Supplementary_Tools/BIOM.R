
# Date : June 15, 2016
# Created by: Dayanara Lebron
# Input is all *.species concatenated together with first column being the sample name
# Purpose: This script creates the BIOM Tables ( OTU_Reads, OTU_RA, and Tax_Table) for input in Phyloseq
# OTU_Reads is a matrix that contains the number of reads per taxid and sample,
# OTU_RA is the same as the OTU_Reads but with relative abundance based on the sample,
# TAX_Table is a matrix that contains all the lineage information for a specific taxid

###Input info from command line####
#! /bin/Rscript
args <- commandArgs(TRUE)
print(args)
filename <- args[1];

#######################Methods################
create_BIOM <- function(x,coltax, colID, colRC) {
  sampleID <- levels(as.factor(x[,colID]))
  
  NCBI_ID <- x[!duplicated(x[,coltax]),coltax]
  
  OTU_TR <- matrix(nrow = length(NCBI_ID), ncol = length(sampleID))
  colnames(OTU_TR) <- sampleID;
  rownames(OTU_TR) <- NCBI_ID;
  
  i <- 1
  while (i <= dim(x)[1])
  {
    for (j in 1:length(sampleID))
    {
      if (identical(x[i,colID],sampleID[j]) == TRUE) {
        break;
      }
    }
    for (k in 1:length(NCBI_ID))
    {
      if (identical(x[i,coltax],NCBI_ID[k]) == TRUE)
      {
        OTU_TR[k,j] <- x[i,colRC]; ##Specify column of Read Counts
      }
    }
    i <- i + 1
  }
  
  for (i in 1:dim(OTU_TR)[1])
  {
    for (j in 1:dim(OTU_TR)[2])
    {
      if (is.na(OTU_TR[i,j]))
        OTU_TR[i,j] <- as.numeric(0)
    }
  }
  
  return(OTU_TR);
}

#Creates Taxonomic Matrix only for Genus and Species

create_TAX <- function(x, coltax, colSpName) {
  Species <- cbind(x[,coltax], x[,colSpName])
  Species <- Species[!duplicated(Species[,1]),]
  rownames(Species) <- Species[,1]
  colnames(Species) <- c("genus","specie")
  
  char <- c('[',']')
  for (i in 1:length(Species[,2]))
  {
    Name <- Species[i,2]
    
    if (grepl("species,", Name))
    {
      tmp <- nchar(Name)
      Species[i,2] <- substr(Name, start = 9, stop = tmp)
    }
    for (n in 1:length(char))
    {
      Species[i,2] <- sub(char[n],"",Species[i,2], fixed = T)
    }
    
  }
  
  for (i in 1:dim(Species)[1])
  {
    Species[i,1] <- strsplit(Species[i,2]," ")[[1]][1]
  }
  return(Species)
}

#Create Contamination Matrix #Adviced only for Human Reads

create_CONT <- function() {
  C_lst <-
    read.table(
      "/usr/gapps/kpath/bin/CONT_TAXID.ls", head = F,stringsAsFactors = F
    )
  OTU_Cont <- OTU_Reads[which(rownames(OTU_Reads) %in% C_lst$V1),]
  TAX_Cont <- Tax_Full[which(rownames(Tax_Full) %in% C_lst$V1),]
}

######################################################################################

spfile <- read.csv(filename, head = F, sep = "\t", stringsAsFactors = F)

rm <- c(9606,32630) #Remove Human and Synthetic construct
spfile <- spfile[-which(spfile[,6] %in% rm),]

#Threshold

threshold = 0;
if (args[2] != 0) {
  threshold = args[2]
}else{
  threshold = 0
}
spfile <- spfile[which(spfile[,3] >= threshold),]

#Create all tables
OTU_Reads <-
  create_BIOM(spfile,6,1,5);  Specie <- create_TAX(spfile,6,7);
OTU_RA <- OTU_Reads;
for (i in 1:dim(OTU_Reads)[2]) {
  OTU_RA[,i] <- OTU_Reads[,i] / sum(OTU_Reads[,i])
}


#Creates Taxonomic Matrix entirely based on the ncbi_taxfile for the current lmat version 
filetax = "./Tax_Ref"
tax <- read.csv(filetax,header = T, sep = "\t")
tax <- tax[which(tax$TAXID %in% rownames(OTU_Reads)),]

if (any(duplicated(tax$TAXID)) == TRUE) {
  tax <- tax[which(!duplicated(tax$TAXID)),]
}

rownames(tax) <- tax$TAXID
Tax_Full <- tax[match(rownames(OTU_Reads),rownames(tax)),]
Tax_Full <- Tax_Full[,-1]


#Identify Contaminants in Sample
#create_Cont();
save(OTU_Reads,Tax_Full,OTU_RA,Specie, file = paste("BIOM_",filename,".RData", sep =""))
