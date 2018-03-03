# Livermore  Metagenomics  Analysis  Toolkit
**Taxonomic classification, content summarization
 and gene identification, an all-in-1 metagenomic analysis toolkit**
___
### Overview

LMAT's main goal is to efficiently assign taxonomic labels to the reads with reference representation down to the species level while maintaining accuracy in the presence of novel organisms. Scalable performance is demonstrated on real and simulated data to show accurate classification even with novel genomes on samples that include viruses, prokaryotes, fungi and protists.

LMAT has three related subcomponents (taxonomic profiling, content summarization and gene annotation) that can be run separately.

### Quick installation

The quick installation procedure will use CMake to ease the process, by downloading, building and installing all the required packages.

#### Required software

* CMake3
* C/C++ compiler with OpenMP support (like gcc, clang, icc, xlc)
* Recommended: python, for some tools
* Optional: MPI, for use in building a Reference Database

#### Using redoall to build LMAT easily

``redoall`` is a convenient wrapper that will direct the installation through CMake for typical compilers (gcc, clang, Intel C/C++ compilers and IBM XL Compilers for Power9):
```
usage: redoall [profile] [compiler]
```

The 1st optional parameter chooses the build profile of CMake:
* ``D`` for ``Debug``
* ``R`` for ``Release`` (this is the current default)
* ``I`` for release with debug info (``RelWithDebInfo``)
* ``M`` for release with minimum size (``MinSizeRel``)
* for just cleaning the parameter is ``clean``

The 2nd optional parameter selects the compiler family:
*  ``gnu`` for using GCC
*  ``intel`` for using Intel compilers
*  ``clang`` for using clang compilers
*  ``ibmpwr9`` for compiling in Power 9 with IBM compilers

#### Example for gcc (release profile)

```
git clone https://github.com/LivGen/LMAT.git
cd LMAT
./redoall
```

### Details

* Scientific details about LMAT are explained in the article "[Scalable metagenomic taxonomy classification using a reference genome database](https://doi.org/10.1093/bioinformatics/btt389)".
* Further information about LMAT can be found in the article "[Using populations of human and microbial genomes for organism detection in metagenomes](https://doi.org/10.1101/gr.184879.114)".
* Please refer to [documentation in the 'doc' subdirectory
  for technical information on LMAT](https://rawgit.com/LivGen/LMAT/master/doc/lmat-doc/index.html).
* [LMAT web site at LLNL](https://computation.llnl.gov/projects/livermore-metagenomics-analysis-toolkit).
* This is an [example of LMAT run](https://sourceforge.net/p/lmat/wiki/Example%20LMAT%20Run/).

### Post-processing with Recentrifuge

If you are analyzing more than one sample with LMAT you can easily visualize and compare them using [Recentrifuge: Robust comparative analysis and contamination removal for metagenomic data](https://github.com/khyox/recentrifuge).

With a score-oriented approach, Recentrifuge is especially useful in the case of low biomass metagenomic studies or when a more reliable detection of minority organisms is needed, like in clinical, environmental and forensic analysis. Further details are in the [bioRxiv pre-print](https://doi.org/10.1101/190934).

For usage and documentation, please, see [running Recentrifuge for LMAT](https://github.com/khyox/recentrifuge/wiki/Running-recentrifuge-for-LMAT) in the [Recentrifuge wiki](https://github.com/khyox/recentrifuge/wiki/). 

### LMAT and PERM: tunning the kernel

LMAT uses [PERM](https://computation.llnl.gov/projects/memory-centric-architectures/perm), a ‘C’ library for persistent heap management used with a dynamic-memory allocator, also developed at LLNL. To PERM (so LMAT) to work in the right conditions, some kernel tunning is advisable:

* Turn off periodic flush to file and dirty ratio flush:
```
echo 0 > /proc/sys/vm/dirty_writeback_centisecs
echo 100 > /proc/sys/vm/dirty_background_ratio
echo 100 > /proc/sys/vm/dirty_ratio
```
* Turn off address space randomization:
```
echo 0 > /proc/sys/kernel/randomize_va_space 
```
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
