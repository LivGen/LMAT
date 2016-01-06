#!/bin/bash
# Proper header for a Bash script.

#
# Rebuild all (or cleans) in LMAT for GNU
#
# You can provide a parameter for choosing the build profile of CMake:
# 	D for Debug (default!)
#	R for Release
#	I for RelWithDebInfo
#	M for MinSizeRel
# For just cleaning the parameter is 'clean'
#
# JMMM - July 2015 - Rel. 0.1
#

# This is specific to LC
echo "Preparing environment in LC clusters to GNU..."
use gcc-4.9.2p
#module -s load gcc/gcc-4.9.0 cmake/3.0.0 openmpi/1.8.1/gnu 
#module load -s python/3.4.0 multinest/3.7/gnu gsl/1.16/gnu
#module list

# This is general
comp="-D COMPILER_FAMILY=gnu" 
echo "Performing the total reconf and rebuild..."
make clean
rm -Rf *.cmake Makefile src/Makefile src/CMakeCache.txt src/*.cmake CMakeCache.txt src/kmerdb/all_headers.hpp CMakeFiles/ src/CMakeFiles/ third-party/*
if [ -n "$1" ]
then
  if [ $1 = "clean" ]
  then
    echo "  Just cleaning, nothing built!"
    exit 0
  elif [ $1 = "R" ]
  then
    echo "  CMake profile will be Release"
    cmake $comp -D CMAKE_BUILD_TYPE=Release .
  elif [ $1 = "I" ]
  then
    echo "  CMake profile will be RelWithDebInfo"
    cmake $comp -D CMAKE_BUILD_TYPE=RelWithDebInfo .
  elif [ $1 = "M" ]
  then
    echo "  CMake profile will be MinSizeRel"
    cmake $comp -D CMAKE_BUILD_TYPE=MinSizeRel .
  elif [ $1 = "D" ]
  then
    echo "  CMake profile will be Debug"
    cmake $comp -D CMAKE_BUILD_TYPE=Debug .
  fi
else
  echo "  CMake profile will be the default one"
  cmake $comp .
fi
make 
