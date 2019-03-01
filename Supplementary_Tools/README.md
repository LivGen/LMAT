## BIOM.R

This script/module converts LMAT output in a BIOM* format that is used in the phyloseq 
R package.  Phyloseq is an R package with nice graphical features to analyze microbial census data.

Some examples can be found here:

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217

and the package can be downloaded here: 

https://bioconductor.org/packages/release/bioc/html/phyloseq.html

### What it does?
The BIOM tables are the following: OTU_Reads, OTU_RA, Tax_Full
-- Lets say that we have a batch of 8 samples, each sampled is processed through LMAT for metagenomic classification. Each sample has n1,...,n8 taxonomical hits, where all $n_{j}$ may not be equal. Then to create an even matrix among all 8 samples the dimension of the OTU_READS=(max{$n_j},8). For each taxonomy in $n_{i}$ not in $n_{j}$ then the Reads are input as 0.



### Run Settings



RScript BIOM.R concatenated.file 

