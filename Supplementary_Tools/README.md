## BIOM.R

This script/module converts LMAT output in a BIOM* format that is used in the phyloseq 
R package.  Phyloseq is an R package with nice graphical features to analyze microbial census data.

Some examples can be found here:

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217

and the package can be downloaded here: 

https://bioconductor.org/packages/release/bioc/html/phyloseq.html

### What it does?
The BIOM tables are the following: OTU_Reads, OTU_RA, Tax_Full
-- Lets say that we have a batch of 8 samples, each sampled is processed through LMAT for metagenomic classification. Each sample has n1,...,n8 taxonomical hits, where all $$n_{j}$$ may not be equal. Then to create an even matrix among all 8 samples the dimension of the OTU_READS=($ max{$n_j} $,8). For each taxonomy in $$n_{i}$$ not in $ n_{j} $ then the Reads are input as 0.

Each column would represent a sample and each row a taxonomy. The OTU_RA works the same but elements are relative abundance of each taxonomical hit per sample. The Tax_Full matrix provides the lineage for each taxid, which is useful to represent the sample per family, genus, etc; depending on whats specified in the phyloseq method.


### Run Settings
Input:
- **concatenated.file** = is all of the *.species* summary outputs per sample of LMAT concatenated together in a single file such that each line has the sample where they belong to.

- **Tax_Ref**= this is a mapping file to obtain the lineage information for a specific taxid based on the reference genomes that are stored in the current LMAT database. [This file is too big to upload, but the script used create it will be made public]

#### Command Line 
RScript BIOM.R $concatenated.file 

Output: 
All the files are stored in a RData environment as: $concatenated.file.RData

Some examples can be seen in the **Example** folder