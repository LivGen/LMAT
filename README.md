# Livermore  Metagenomics  Analysis  Toolkit
**Taxonomic scoring and classification, content summarization
 and gene identification, an all-in-one metagenomic analysis toolkit**
___
### Overview

LMAT's main goal is to efficiently assign taxonomic labels to the reads with reference representation down to the species level but maintain accuracy in the presence of novel organisms.

Scalable performance is demonstrated on real and simulated data to show accurate classification even with novel genomes on samples that include viruses, prokaryotes, fungi and protists.

Each part of LMAT (taxonomic profiling, content summarization and gene annotation) can be run independently.

### Quick installation

The quick installation will use CMake to ease the process, by downloading, building and installing all the required packages.

#### Required software

* CMake3
* C/C++ compiler with OpenMP support (like gcc, clang, icc)
* Recommended: python, for some tools
* Optional: MPI, for use in building a Reference Database

#### Using the wrappers around CMake

There are some wrappers that will direct the installation through CMake for typical compilers (gcc, clang and Intel C/C++ compilers). These wrapers accept a single parameter for the type of build:
* ``D`` for ``Debug`` (this is the default one)
* ``R`` for ``Release``
* ``I`` for release with debug info (``RelWithDebInfo``)
* ``M`` for release with minimum size (``MinSizeRel``)
* for just cleaning the parameter is ``clean``

#### Example for gcc

```
git clone https://github.com/LivGen/LMAT.git
./redoall_gnu.sh
```

### Details

* Scientific details about LMAT are explained in the article "[Scalable metagenomic taxonomy classification using a reference genome database](https://doi.org/10.1093/bioinformatics/btt389)".
* Further information about LMAT can be found in the article "[Using populations of human and microbial genomes for organism detection in metagenomes](https://doi.org/10.1101/gr.184879.114)".
* Please refer to [documentation in the 'doc' subdirectory
  for technical information on LMAT](https://rawgit.com/LivGen/LMAT/master/doc/lmat-doc/index.html).


### Post-processing with Recentrifuge

If you are analyzing more than one sample with LMAT you can easily analyze all of them using [Recentrifuge: Robust comparative analysis and contamination removal for metagenomic data](https://github.com/khyox/recentrifuge).

Recentrifuge is especially useful in the case of low biomass metagenomic studies or when a more reliable detection of minority organisms is needed, like in clinical, environmental and forensic analysis. Further details in the [bioRxiv pre-print](https://doi.org/10.1101/190934).

For usage and documentation, please, see [the Recentrifuge wiki](https://github.com/khyox/recentrifuge/wiki).

___

```
 ==============================================
  : : :       ··        ··      ·   ·········· 
  : : :       ···      ···     · ·      ··     
  : : :       ·· ··  ·· ··    ·· ··     ··     
  : : ······· ··   ··   ··   ·······    ··     
  :  ······   ··        ··  ··     ··   ··     
    ·····     ··        ·· ··       ··  ··     
 ==============================================
   Livermore  Metagenomics  Analysis  Toolkit  
 ==============================================
```




